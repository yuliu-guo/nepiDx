import os
import argparse
import pandas as pd
import subprocess
import re
import io
import csv

def load_p_group_map(p_group_file):
    """
    Parses the P group file and returns a mapping of full allele (e.g. 01:315)
    to its P group root (e.g. 01:03), using gene prefix (e.g. A*).
    Skips entries with missing P-group root values.
    """
    allele_to_pgroup = {}

    with open(p_group_file, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split(';')
            if len(fields) < 3:
                continue

            gene = fields[0].strip()           # e.g., 'A*'
            alleles = fields[1].split('/')     # e.g., ['01:315', '01:01:153']
            p_group = fields[2].strip()        # e.g., '01:03P'

            if not p_group:                    # Skip if no P group root
                continue

            # Extract root from p_group (removing trailing 'P')
            root = p_group.replace('P', '')

            for a in alleles:
                trimmed = a.strip()
                if trimmed:
                    allele_to_pgroup[f'{gene}{trimmed}'] = f'{gene}{root}'

    return allele_to_pgroup


def map_to_p_group(full_allele, pgroup_map):
    """
    Maps full allele to its HLA P-group root (e.g., HLA-A*01:315 -> HLA-A*01:03),
    or returns the original allele if no mapping exists.
    """
    match = re.match(r'^(HLA-[A-Z]\*)(\S+)', full_allele)
    if not match:
        return full_allele  # fallback - return full allele even if pattern is weird

    prefix, allele_part = match.groups()       # e.g., 'HLA-A*', '01:12'
    pgroup_key = f'{prefix[4:]}{allele_part}'  # strip 'HLA-' -> 'A*01:12'
    
    if pgroup_key in pgroup_map:
        return f'HLA-{pgroup_map[pgroup_key]}'
    else:
        return full_allele                      # no mapping found, return as-is

def extract_and_trim_alleles(hla_allele_path, p_group_file, hla_quality_threshold=0):
    """
    Reads allele data and returns list of trimmed alleles with quality >= threshold.
    """
    hla_pattern = re.compile(r'^(HLA-[ABC]\*\S+)\s+(\d+)$')
    processed_alleles = []
    trimmed_alleles = []

    # Load P group mapping
    pgroup_map = load_p_group_map(p_group_file)

    with open(hla_allele_path, 'r') as file:
        for line in file:
            match = hla_pattern.match(line.strip())
            if match:
                allele, quality = match.group(1), int(match.group(2))
                if quality >= hla_quality_threshold:
                    trimmed = map_to_p_group(allele, pgroup_map)
                    if trimmed:
                        trimmed_alleles.append(trimmed)
                        processed = trimmed.replace(':', '').replace('*', '')
                        processed_alleles.append(processed)
    alleles_df = pd.DataFrame({
        "Full_Allele": trimmed_alleles,
        "Processed_Allele": processed_alleles
    })
    return alleles_df

def get_available_alleles(allele_list_file):
    """
    Reads the available alleles from a file and returns a list of allele names.
    """
    available_alleles = []
    with open(allele_list_file, 'r') as f:
        for line in f:
            available_alleles.append(line.strip())
    return available_alleles

def run_ICERFIRE_prediction(patient_id, input_peptides, input_alleles, available_alleles_file, output_file, scratch_dir, p_group_file, hla_quality_threshold=0):
    """
    Runs ICERFIRE prediction with given peptides and HLA alleles and saves the output to a .csv file.
    Parameters:
        patient_id (str): The ID of the patient for whom to run predictions
        input_peptides (str): A path to the input file containing peptides
        input_alleles (str): A path to the input file containing HLA alleles
        output_file (str): A path to the output .csv file
        scratch_dir (str): A path to the scratch directory for temporary files
    """
    ### Read peptides from similarity file and create input for ICERFIRE
    peptides_df = pd.read_csv(input_peptides, sep="\t")
    peptides_df = peptides_df[["ID", "Mut_peptide", "WT_peptide"]]
    peptides_df = peptides_df[peptides_df['Mut_peptide'] != peptides_df['WT_peptide']]
    
    valid_amino_acids = r'^[ACDEFGHIKLMNPQRSTVWY]+$'
    
    is_valid_mut = peptides_df['Mut_peptide'].str.match(valid_amino_acids, na=False)
    is_valid_wt = peptides_df['WT_peptide'].str.match(valid_amino_acids, na=False)
    
    peptides_df = peptides_df[is_valid_mut & is_valid_wt]
    
    ### STOP if there are no peptides
    if len(peptides_df) == 0:
        print("No peptides available in the nmer file.")
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow([
                  'Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Immunogenicity_Score', 'Immunogenicity_Rank', 'Immunogenicity_tool'
            ]) #write an empty header row
        return 
    
    ### Read HLA alleles
    alleles_df = extract_and_trim_alleles(input_alleles, p_group_file, hla_quality_threshold)
    
    ### Get available alleles for ICERFIRE prediction
    all_available_alleles = get_available_alleles(available_alleles_file)
    
    skipped_alleles = alleles_df[~alleles_df['Processed_Allele'].isin(all_available_alleles)]['Full_Allele'].tolist()
    
    ### STOP if there are no alleles
    if len(alleles_df[alleles_df['Processed_Allele'].isin(all_available_alleles)]) == 0:
        print("No supported alleles passed the quality threshold.")
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow([
                  'Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Immunogenicity_Score', 'Immunogenicity_tool'
            ]) #write an empty header row
        return
    
    peptides_alleles_df = pd.DataFrame([(peptide_id, peptide, wt_peptide, allele_processed, allele_full) 
                                        for peptide_id, peptide, wt_peptide in zip(peptides_df['ID'], peptides_df['Mut_peptide'], peptides_df['WT_peptide']) 
                                        for allele_processed, allele_full in zip(alleles_df['Processed_Allele'], alleles_df['Full_Allele'])], 
                                       columns=['Identity', 'Peptide', 'WT_Peptide', 'Processed_Allele', 'Full_Allele'])
    peptides_alleles_filtered_df = peptides_alleles_df[~peptides_alleles_df['Full_Allele'].isin(skipped_alleles)]
    
    ICERFIRE_input = peptides_alleles_filtered_df[['Peptide', 'WT_Peptide', 'Processed_Allele']].copy().drop_duplicates()
    ICERFIRE_input.to_csv(f"{scratch_dir}/ICERFIRE_input.txt", sep=",", index=False, header=False)
    
    ## Run ICERFIRE
    command = (
        f'ICERFIRE.sh -f {scratch_dir}/ICERFIRE_input.txt -a false -u false -o {scratch_dir}/ICERFIRE_tmp_output.csv'
    )
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    
    ### Print the STDOUT and STDERR for debugginf purposes
    print(result.stdout)
    print(result.stderr)
    result_df = pd.read_csv(f"{scratch_dir}/ICERFIRE_tmp_output.csv")
    result_df = result_df.rename(columns={'HLA': 'Processed_Allele', 
                                          'wild_type': 'WT_Peptide', 
                                          'prediction': 'Immunogenicity_Score',
                                          '%Rank': 'Immunogenicity_Rank'})
    result_df['Processed_Allele'] = result_df['Processed_Allele'].str.replace('*', '')
    
    result_df['Immunogenicity_tool'] = 'ICERFIRE'
    result_df['Patient_ID'] = patient_id
    result_df = result_df[['Patient_ID', 'Processed_Allele', 'Peptide', 'WT_Peptide', 'Immunogenicity_Score', 'Immunogenicity_Rank', 'Immunogenicity_tool']]
    result_df = pd.merge(peptides_alleles_filtered_df, result_df, on=['Peptide', 'WT_Peptide', 'Processed_Allele'], how='left')
    result_df = result_df.rename(columns={'Full_Allele': 'HLA_allele'})
    result_df = result_df[['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Immunogenicity_Score', 'Immunogenicity_Rank', 'Immunogenicity_tool']]
    
    ### Pairs skipped because allele is not available
    peptides_alleles_skipped_df = peptides_alleles_df[peptides_alleles_df['Full_Allele'].isin(skipped_alleles)]
    
    peptides_alleles_skipped_df['Patient_ID'] = patient_id
    peptides_alleles_skipped_df['Immunogenicity_Score'] = None
    peptides_alleles_skipped_df['Immunogenicity_Rank'] = None
    peptides_alleles_skipped_df['Immunogenicity_tool'] = 'ICERFIRE'
    peptides_alleles_skipped_df = peptides_alleles_skipped_df.rename(columns={'Full_Allele': 'HLA_allele'})
    peptides_alleles_skipped_df = peptides_alleles_skipped_df[['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Immunogenicity_Score', 'Immunogenicity_Rank', 'Immunogenicity_tool']]
    
    result_df = pd.concat([result_df, peptides_alleles_skipped_df], ignore_index=True)
    result_df.drop_duplicates(inplace=True)
    result_df.to_csv(output_file, index=False, sep="\t")
    

def cleanup(scratch_dir):
    """
    Cleans up temporary files in the scratch directory.
    """
    temp_input = os.path.join(scratch_dir, "ICERFIRE_input.txt")
    temp_output = os.path.join(scratch_dir, "ICERFIRE_tmp_output.csv")
    if os.path.exists(temp_output):
        os.remove(temp_output)
    if os.path.exists(temp_input):
        os.remove(temp_input)

def main():
    parser = argparse.ArgumentParser(description="Run ICERFIRE predictions")
    parser.add_argument("--patient_id", required=True, help="The patient ID for the predictions")
    parser.add_argument("--input_peptides", required=True, help="Path to the input nmer peptides file (.fasta)")
    parser.add_argument("--input_alleles", required=True, help="Path to the input nmer alleles file (.fasta)")
    parser.add_argument("--alleles_list", required=True, help="Path to the list of alleles available for IEDB prediction.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to save the output file")
    parser.add_argument("--scratch_dir", type=str, required=True, help="Path to the scratch directory")
    parser.add_argument("--p_group_file", type=str, required=True, help="Path to P_GROUP_FILE for HLA allele mapping")
    parser.add_argument("--hla_quality_threshold", type=int, required=True, help="Quality threshold for filtering alleles")

    args = parser.parse_args()

    try:
        run_ICERFIRE_prediction(args.patient_id, args.input_peptides, args.input_alleles, args.alleles_list, args.output_file, args.scratch_dir, args.p_group_file, args.hla_quality_threshold)
    finally:
        #cleanup(args.scratch_dir)
        pass

if __name__ == "__main__":
    main()