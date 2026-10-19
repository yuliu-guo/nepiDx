import os
import argparse
import pandas as pd
import subprocess
import re
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
        final_allele = full_allele  # fallback - full allele even if pattern is weird
    else:
        prefix, allele_part = match.groups()       # e.g., 'HLA-A*', '01:12'
        pgroup_key = f'{prefix[4:]}{allele_part}'  # strip 'HLA-' -> 'A*01:12'
        
        if pgroup_key in pgroup_map:
            final_allele = f'HLA-{pgroup_map[pgroup_key]}'
        else:
            final_allele = full_allele                      # no mapping found
            
    # Validity Check
    # An allele is valid only if the part after the '*' contains no alphabetical letters
    if '*' in final_allele:
        allele_designation = final_allele.split('*', 1)[1]
        # If any letter (a-z, A-Z) is found in the numbers/digits group, it's invalid
        valid = not bool(re.search(r'[a-zA-Z]', allele_designation))
    else:
        # If the format doesn't even contain a '*', flag it as invalid
        valid = False

    return final_allele, valid

def read_fasta_to_df(filepath):
    """
    Reads a FASTA file and returns a DataFrame with columns 'Identity' and 'Peptide'.
    Assumes that each FASTA sequence is in one line.
    """
    
    with open(filepath, 'r') as f:
        # Read lines and strip whitespace/newlines
        lines = [line.strip() for line in f if line.strip()]

    data = {
        "Identity": [l[1:] for l in lines[0::2]], 
        "Peptide": lines[1::2]
    }
    
    return pd.DataFrame(data)

def extract_and_trim_alleles(hla_allele_path, p_group_file, hla_quality_threshold=0):
    """
    Reads allele data and returns list of trimmed alleles with quality >= threshold.
    """
    hla_pattern = re.compile(r'^(HLA-[ABC]\*\S+)\s+(\d+)$')
    processed_alleles = []
    trimmed_alleles = []
    is_allele_valid = []

    # Load P group mapping
    pgroup_map = load_p_group_map(p_group_file)

    with open(hla_allele_path, 'r') as file:
        for line in file:
            match = hla_pattern.match(line.strip())
            if match:
                allele, quality = match.group(1), int(match.group(2))
                if quality >= hla_quality_threshold:
                    trimmed, valid = map_to_p_group(allele, pgroup_map)
                    if trimmed:
                        trimmed_alleles.append(trimmed)
                        processed = trimmed.replace(':', '')
                        processed_alleles.append(processed)
                        is_allele_valid.append(valid)
                        
    alleles_df = pd.DataFrame({
        "Full_Allele": trimmed_alleles,
        "Processed_Allele": processed_alleles,
        "Is_valid": is_allele_valid
    })
    return alleles_df

def run_deepimmuno_prediction(patient_id, input_peptides, input_alleles, output_file, scratch_dir, p_group_file, hla_quality_threshold=0):
    """
    Runs DeepImmuno prediction with given peptides and HLA alleles and saves the output to a .csv file.
    Parameters:
        patient_id (str): The patient ID for the predictions
        input_peptides (str): A path to the input file containing peptides
        input_alleles (str): A path to the input file containing HLA alleles
        output_file (str): A path to the output .csv file
        scratch_dir (str): A path to the scratch directory for temporary files
    """
    ### Read peptides from FASTA file
    peptides_df = read_fasta_to_df(input_peptides)
    
    valid_amino_acids = r'^[ACDEFGHIKLMNPQRSTVWY]+$'
    
    is_valid = peptides_df['Peptide'].str.match(valid_amino_acids, na=False)
    peptides_df = peptides_df[is_valid]
    
    ### STOP if there are no peptides
    if len(peptides_df) == 0:
        print("No peptides available in the nmer file.")
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow([
                  'Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Immunogenicity_Score', 'Immunogenicity_tool'
            ]) #write an empty header row
        return 
    
    ### Read HLA alleles
    alleles_df = extract_and_trim_alleles(input_alleles, p_group_file, hla_quality_threshold)
    
    ### STOP if there are no alleles
    if len(alleles_df[alleles_df['Is_valid']]) == 0:
        print("No valid alleles passed the quality threshold.")
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow([
                  'Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Immunogenicity_Score', 'Immunogenicity_tool'
            ]) #write an empty header row
        return 
    
    skipped_alleles = alleles_df[~alleles_df['Is_valid']]['Full_Allele'].tolist()
    
    ### Prepare input for DeepImmuno
    peptides_alleles_df = pd.DataFrame([(peptide_id, peptide, allele_processed, allele_full) 
                                        for peptide_id, peptide in zip(peptides_df['Identity'], peptides_df['Peptide']) 
                                        for allele_processed, allele_full in zip(alleles_df['Processed_Allele'], alleles_df['Full_Allele'])], 
                                       columns=['Identity', 'Peptide', 'Processed_Allele', 'Full_Allele'])
                                       
    peptides_alleles_filtered_df = peptides_alleles_df[~peptides_alleles_df['Full_Allele'].isin(skipped_alleles)]
    
    deepimmuno_input = peptides_alleles_filtered_df[['Peptide', 'Processed_Allele']]
    deepimmuno_input.to_csv(f"{scratch_dir}/deepimmuno_input.csv", sep=",", index=False, header=False)

    ## Run DeepImmuno
    command = (
        f'python3 /deepimmuno-cnn.py --mode "multiple" --intdir "{scratch_dir}/deepimmuno_input.csv" --outdir "{scratch_dir}"'
    )
    try:
        result = subprocess.run(command, check=True, capture_output=False, text=True, shell=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to run DeepImmuno: {e}")
    
    ## Process DeepImmuno output
    deepimmuno_df = pd.read_csv(f"{scratch_dir}/deepimmuno-cnn-result.txt", sep="\t")
    deepimmuno_df = deepimmuno_df.rename(columns={'peptide': 'Peptide', 'HLA': 'Processed_Allele', 'immunogenicity': 'Immunogenicity_Score'})
    deepimmuno_df = deepimmuno_df.drop_duplicates(subset=['Peptide', 'Processed_Allele'])
    result_df = pd.merge(peptides_alleles_filtered_df, deepimmuno_df, on=['Peptide', 'Processed_Allele'], how='left')
    result_df['Immunogenicity_tool'] = 'DeepImmuno'
    result_df = result_df.drop(columns=['Processed_Allele'])
    result_df = result_df.rename(columns={'Full_Allele': 'HLA_allele'})
    result_df.insert(0, 'Patient_ID', patient_id)
    result_df = result_df[['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Immunogenicity_Score', 'Immunogenicity_tool']]
    
    ### Pairs skipped because allele is not valid
    peptides_alleles_skipped_df = peptides_alleles_df[peptides_alleles_df['Full_Allele'].isin(skipped_alleles)].copy()
    peptides_alleles_skipped_df['Patient_ID'] = patient_id
    peptides_alleles_skipped_df['Immunogenicity_Score'] = None
    peptides_alleles_skipped_df['Immunogenicity_tool'] = 'DeepImmuno'
    peptides_alleles_skipped_df = peptides_alleles_skipped_df.rename(columns={'Full_Allele': 'HLA_allele'})
    peptides_alleles_skipped_df = peptides_alleles_skipped_df[['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Immunogenicity_Score', 'Immunogenicity_tool']]
    
    result_df = pd.concat([result_df, peptides_alleles_skipped_df], ignore_index=True)
    result_df.to_csv(output_file, sep="\t", index=False)

def cleanup(scratch_dir):
    """
    Cleans up temporary files in the scratch directory.
    """
    temp_input = os.path.join(scratch_dir, "deepimmuno_input.csv")
    temp_result = os.path.join(scratch_dir, "deepimmuno-cnn-result.txt")
    if os.path.exists(temp_input):
        os.remove(temp_input)
    if os.path.exists(temp_result):
        os.remove(temp_result)

def main():
    parser = argparse.ArgumentParser(description="Run DeepImmuno predictions")
    parser.add_argument("--patient_id", required=True, help="The patient ID for the predictions")
    parser.add_argument("--input_peptides", required=True, help="Path to the input nmer peptides file (.fasta)")
    parser.add_argument("--input_alleles", required=True, help="Path to the input nmer alleles file (.fasta)")
    parser.add_argument("-o", "--output_file", required=True, help="Path to save the output file")
    parser.add_argument("--scratch_dir", type=str, required=True, help="Path to the scratch directory")
    parser.add_argument("--p_group_file", type=str, required=True, help="Path to P_GROUP_FILE for HLA allele mapping")
    parser.add_argument("--hla_quality_threshold", type=int, required=True, help="Quality threshold for filtering alleles")

    args = parser.parse_args()

    try:
        run_deepimmuno_prediction(args.patient_id, args.input_peptides, args.input_alleles, args.output_file, args.scratch_dir, args.p_group_file, args.hla_quality_threshold)
    finally:
        cleanup(args.scratch_dir)

if __name__ == "__main__":
    main()