import os
import subprocess
import argparse
import pandas as pd
import csv
import re

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
        return full_allele  # fallback — return full allele even if pattern is weird

    prefix, allele_part = match.groups()       # e.g., 'HLA-A*', '01:12'
    pgroup_key = f'{prefix[4:]}{allele_part}'  # strip 'HLA-' -> 'A*01:12'
    
    if pgroup_key in pgroup_map:
        return f'HLA-{pgroup_map[pgroup_key]}'
    else:
        return full_allele                      # no mapping found, return as-is
        
def get_first_peptide_length(fasta_path):
    """
    Reads a FASTA file and returns the length of the first sequence.
    Parameters:
        fasta_path (str): A path to the FASTA file
    Returns:
        int: Length of the first peptide sequence in the file
    """
    with open(fasta_path) as f:
        seq_lines = []
        in_first_record = False
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if in_first_record:
                    # we've already collected the first record's sequence; stop
                    break
                in_first_record = True
                continue
            if in_first_record:
                seq_lines.append(line)

    if not seq_lines:
        raise ValueError(f"No sequence found in FASTA file: {fasta_path}")

    return len(''.join(seq_lines))

def extract_and_trim_alleles(hla_allele_path, p_group_file, hla_quality_threshold=0):
    """
    Reads allele data and returns list of trimmed alleles with quality >= threshold.
    """
    hla_pattern = re.compile(r'^(HLA-[ABC]\*\S+)\s+(\d+)$')
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

    return trimmed_alleles

def run_mhcflurry(patient_id, trimmed_alleles, mutated_fasta_path, output_file, scratch_dir):
    """
    Runs MHCflurry prediction with given alleles and peptides and saves the output to a .csv file.

    Parameters:
        patient_id (str): The patient identifier
        trimmed_alleles (list): A list of trimmed HLA alleles
        mutated_fasta_path (str): A path to the mutated peptides FASTA file
        output_file (str): A path to the output .csv file
        scratch_dir (str): A path to a directory for storing temporary files; the files are removed after the process finishes
    """
    peptide_length = get_first_peptide_length(mutated_fasta_path)
    
    df_results = pd.DataFrame(columns=['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Affinity', 'Affinity_percentile', 'Processing_score', 'Presentation_score', 'Presentation_percentile', 'Prediction_tool'])
    for allele in trimmed_alleles:
      command = (
        f"mhcflurry-predict-scan {mutated_fasta_path} --alleles {allele} --out {scratch_dir}/mhcflurry_output.csv --results-all --peptide-lengths {peptide_length} --input-format fasta --no-throw"
        )
      try:
        result = subprocess.run(command, check=True, capture_output=True, text=True, shell=True)
        df_tmp = pd.read_csv(f"{scratch_dir}/mhcflurry_output.csv")
        df_tmp = df_tmp[['best_allele', 'peptide', 'sequence_name', 'affinity', 'affinity_percentile', 'processing_score', 'presentation_score', 'presentation_percentile']]
        df_tmp.insert(0, 'Patient_ID', patient_id)
        df_tmp['Prediction_tool'] = 'MHCflurry'
        df_tmp.columns = df_results.columns
        df_results = pd.concat([df_results, df_tmp], ignore_index=True)
      except subprocess.CalledProcessError as e:
        print(f"Failed to run MHCflurry for allele {allele}: {e}")
        continue
      except pandas.errors.EmptyDataError as e:
        print(f"No results for allele {allele}: {e}")
        continue

    df_results.to_csv(output_file, sep = '\t', index=False)

def cleanup(scratch_dir):
    """
    Cleans up temporary files in the scratch directory.
    """
    temp_file = os.path.join(scratch_dir, "mhcflurry_output.csv")
    if os.path.exists(temp_file):
        os.remove(temp_file)

def main():
    parser = argparse.ArgumentParser(description="Run MHCflurry predictions")
    parser.add_argument("-p", "--patient_id", required=True, help="Patient ID")
    parser.add_argument("--hla_allele_path", required=True, help="Path to the file containing HLA alleles and quality scores")
    parser.add_argument("--hla_quality_threshold", type=int, required=True, help="Quality threshold for filtering alleles")
    parser.add_argument("--mutated_peptides_fasta", required=True, help="Path to the mutated peptides FASTA file")
    parser.add_argument("-o", "--output_file", required=True, help="Path to save the output files")
    parser.add_argument("--scratch_dir", type=str, required=True, help="Path to the scratch directory")
    parser.add_argument("--p_group_file", type=str, required=True, help="Path to P_GROUP_FILE for HLA allele mapping")

    args = parser.parse_args()

    try:
        trimmed_alleles = extract_and_trim_alleles(args.hla_allele_path, args.p_group_file, args.hla_quality_threshold)
        if not trimmed_alleles:
            print("No allele passed the quality threshold")
            with open(args.output_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow(['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Affinity', 'Affinity_percentile', 
                            'Processing_score', 'Presentation_score', 'Presentation_percentile', 'Prediction_tool']) #write an empty header row
            return 

        run_mhcflurry(args.patient_id, trimmed_alleles, args.mutated_peptides_fasta, args.output_file, args.scratch_dir)
    finally:
        cleanup(args.scratch_dir)
    

if __name__ == "__main__":
    main()
