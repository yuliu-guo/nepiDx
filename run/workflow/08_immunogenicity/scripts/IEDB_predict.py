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

def run_IEDB_prediction(patient_id, input_peptides, input_alleles, available_alleles_file, output_file, scratch_dir, p_group_file, hla_quality_threshold=0):
    """
    Runs IEDB prediction with given peptides and HLA alleles and saves the output to a .csv file.
    Parameters:
        patient_id (str): The ID of the patient for whom to run predictions
        input_peptides (str): A path to the input file containing peptides
        input_alleles (str): A path to the input file containing HLA alleles
        output_file (str): A path to the output .csv file
        scratch_dir (str): A path to the scratch directory for temporary files
    """
    ### Read peptides from FASTA file and create input for IEDB
    peptides_df = read_fasta_to_df(input_peptides)
    
    valid_amino_acids = r'^[ACDEFGHIKLMNPQRSTVWY]+$'
    
    is_valid = peptides_df['Peptide'].str.match(valid_amino_acids, na=False)
    peptides_df = peptides_df[is_valid]
    
    peptides_list = peptides_df['Peptide'].tolist()
    with open(f"{scratch_dir}/IEDB_peptides.txt", 'w') as f:
        for peptide in peptides_list:
            f.write(f"{peptide}\n")
    
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
    
    ### Get available alleles for IEDB prediction
    available_alleles = get_available_alleles(available_alleles_file)
    
    ### STOP if there are no alleles
    if len(alleles_df[alleles_df['Processed_Allele'].isin(available_alleles)]) == 0:
        print("No supported alleles passed the quality threshold.")
        with open(output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow([
                  'Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Immunogenicity_Score', 'Immunogenicity_tool'
            ]) #write an empty header row
        return
    
    result_df = pd.DataFrame(columns=['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Immunogenicity_Score', 'Immunogenicity_tool'])
    
    ## Run IEDB
    for allele, full_allele in zip(alleles_df['Processed_Allele'], alleles_df['Full_Allele']):
        if allele in available_alleles:
            command = (
                f'python3 /predict_immunogenicity.py "{scratch_dir}/IEDB_peptides.txt" -a {allele}'
            )
            try:
                result = subprocess.run(command, check=True, capture_output=True, text=True, shell=True)
                output = result.stdout
                
                ## Handle IEDB output and save to result_df
                tmp_df = pd.read_csv(io.StringIO(result.stdout), skiprows=4)
                tmp_df = tmp_df[['peptide', 'score']]
                tmp_df = tmp_df.drop_duplicates(subset=['peptide'])
                tmp_df['HLA_allele'] = full_allele
                tmp_df['Immunogenicity_tool'] = 'IEDB'
                tmp_df = tmp_df.rename(columns={'peptide': 'Peptide', 'score': 'Immunogenicity_Score'})
                tmp_df['Patient_ID'] = patient_id
                tmp_df = tmp_df.merge(peptides_df, on='Peptide', how='left')
                tmp_df = tmp_df[['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Immunogenicity_Score', 'Immunogenicity_tool']]
                result_df = pd.concat([result_df, tmp_df], ignore_index=True)

            except subprocess.CalledProcessError as e:
                print(f"Failed to run IEDB: {e}")
        else:
            print(f"Allele {allele} is not available for IEDB prediction.")
            tmp_df = pd.DataFrame({
                'Patient_ID': patient_id,
                'HLA_allele': full_allele,
                'Peptide': peptides_df['Peptide'],
                'Identity': peptides_df['Identity'],
                'Immunogenicity_Score': [None] * len(peptides_df),
                'Immunogenicity_tool': 'IEDB'
            })
            result_df = pd.concat([result_df, tmp_df], ignore_index=True)
            
    result_df.to_csv(output_file, sep="\t", index=False)

def cleanup(scratch_dir):
    """
    Cleans up temporary files in the scratch directory.
    """
    temp_input = os.path.join(scratch_dir, "IEDB_peptides.txt")
    if os.path.exists(temp_input):
        os.remove(temp_input)

def main():
    parser = argparse.ArgumentParser(description="Run IEDB predictions")
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
        run_IEDB_prediction(args.patient_id, args.input_peptides, args.input_alleles, args.alleles_list, args.output_file, args.scratch_dir, args.p_group_file, args.hla_quality_threshold)
    finally:
        cleanup(args.scratch_dir)

if __name__ == "__main__":
    main()