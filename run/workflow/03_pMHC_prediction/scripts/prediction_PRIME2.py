import os
import subprocess
import argparse
import pandas as pd
import csv
import re
from pathlib import Path
from collections import defaultdict

PEPTIDES_FILENAME = "peptides.txt"
PRIME_OUTPUT_FILENAME = "PRIME2_output.csv"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

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
         # mapped value
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
                        processed = trimmed.replace('*', '').replace('-', '').replace(':', '').replace('HLA', '')  # Remove the asterisk
                        processed_alleles.append(processed)

    return processed_alleles, trimmed_alleles

def process_peptides(fasta_path, scratch_dir):
    """
    Reads peptides from a FASTA file and writes all unique peptide sequences
    to a single file in scratch_dir (one per line), as PRIME input.

    Returns:
        seq_to_identities (dict): {sequence: [identity, ...]} lookup, so results for a
                                  sequence can be saved for every identity that has it
    """
    seq_to_identities = defaultdict(list)
    identity, seq_parts = None, []
    
    with open(fasta_path, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if identity is not None:
                    seq_to_identities["".join(seq_parts)].append(identity)
                identity, seq_parts = line[1:].strip(), []
            else:
                seq_parts.append(line)
                
        if identity is not None:
            seq_to_identities["".join(seq_parts)].append(identity)
            
    invalid = [seq for seq in seq_to_identities if set(seq) - VALID_AA]
    for seq in invalid:
        print(f"Skipping peptide {seq} (identities: {seq_to_identities[seq]}): "
              f"non-standard amino acids {sorted(set(seq) - VALID_AA)}")
        del seq_to_identities[seq]
        
    with open(os.path.join(scratch_dir, PEPTIDES_FILENAME), 'w') as f:
        for seq in seq_to_identities:
            f.write(f"{seq}\n")

    return dict(seq_to_identities)

def run_PRIME2(patient_id, processed_alleles, trimmed_alleles, seq_to_identities, output_file, scratch_dir):
    """
    Runs PRIME2 once per allele on all peptides and saves the results to a .csv file,
    with one row per (allele, peptide identity).

    Parameters:
        patient_id (str): The patient identifier
        processed_alleles (list): A list of trimmed and processed HLA alleles
        trimmed_alleles (list): A list of trimmed HLA alleles
        seq_to_identities (dict): {sequence: [identity, ...]} lookup from process_peptides
        output_file (str): A path to the output .csv file
        scratch_dir (str): A path to the directory holding the temporary files
    """
    columns = ['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', '%Rank', 'Score',
               '%RankBinding', 'Prediction_tool']
               
    lookup_df = (
        pd.DataFrame({'Peptide': list(seq_to_identities.keys()),
                      'Identity': list(seq_to_identities.values())})
        .explode('Identity')
    )

    peptides_file = os.path.join(scratch_dir, PEPTIDES_FILENAME)
    prime_output = os.path.join(scratch_dir, PRIME_OUTPUT_FILENAME)
    results = []
    for allele, trimmed_allele in zip(processed_alleles, trimmed_alleles):
        command = (
                f"PRIME -i {peptides_file} -o {prime_output} -a {allele}"
            )
        try:
            subprocess.run(command, check=True, capture_output=True, text=True, shell=True)
            df_tmp = pd.read_csv(prime_output, comment='#', sep='\t')
        except subprocess.CalledProcessError as e:
            print(f"Failed to run PRIME2 for allele {trimmed_allele}: {e}")
            continue

        df_tmp = df_tmp[['Peptide', '%Rank_bestAllele', 'Score_bestAllele',
                         '%RankBinding_bestAllele']]
        df_tmp = df_tmp.merge(lookup_df, on='Peptide', how='inner')
        df_tmp.insert(0, 'Patient_ID', patient_id)
        df_tmp.insert(1, 'HLA_allele', trimmed_allele)
        df_tmp['Prediction_tool'] = 'PRIME2'
        df_tmp = df_tmp[['Patient_ID', 'HLA_allele', 'Peptide', 'Identity',
                         '%Rank_bestAllele', 'Score_bestAllele',
                         '%RankBinding_bestAllele', 'Prediction_tool']]
        df_tmp.columns = columns
        results.append(df_tmp)

    df_results = pd.concat(results, ignore_index=True) if results else pd.DataFrame(columns=columns)
    df_results.to_csv(output_file, sep='\t', index=False)

def cleanup(scratch_dir):
    """
    Cleans up temporary files in the scratch directory.
    """
    for filename in (PEPTIDES_FILENAME, PRIME_OUTPUT_FILENAME):
        path = os.path.join(scratch_dir, filename)
        if os.path.exists(path):
            os.remove(path)

def main():
    parser = argparse.ArgumentParser(description="Run PRIME2 predictions")
    parser.add_argument("-p", "--patient_id", required=True, help="Patient ID")
    parser.add_argument("--hla_allele_path", required=True, help="Path to the file containing HLA alleles and quality scores")
    parser.add_argument("--hla_quality_threshold", type=int, required=True, help="Quality threshold for filtering alleles")
    parser.add_argument("--mutated_peptides_fasta", required=True, help="Path to the mutated peptides FASTA file")
    parser.add_argument("-o", "--output_file", required=True, help="Path to save the output files")
    parser.add_argument("--scratch_dir", type=str, required=True, help="Path to the scratch directory")
    parser.add_argument("--p_group_file", type=str, required=True, help="Path to P_GROUP_FILE for HLA allele mapping")

    args = parser.parse_args()

    try:
        processed_alleles, trimmed_alleles = extract_and_trim_alleles(args.hla_allele_path, args.p_group_file, args.hla_quality_threshold)
        if not trimmed_alleles:
            print("No allele passed the quality threshold")
            with open(args.output_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow(['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', '%Rank', 
                                'Score', '%RankBinding', 'Prediction_tool']) #write an empty header row
            return 

        seq_to_identities = process_peptides(args.mutated_peptides_fasta, args.scratch_dir)
        run_PRIME2(args.patient_id, processed_alleles, trimmed_alleles, 
                    seq_to_identities, args.output_file, args.scratch_dir)
    finally:
        cleanup(args.scratch_dir)
    

if __name__ == "__main__":
    main()