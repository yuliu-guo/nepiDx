import os
import subprocess
import argparse
import pandas as pd
import numpy as np
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
    
def similarity_to_fasta(similarity_file, scratch_dir):
    sim_df = pd.read_csv(similarity_file, sep="\t")
    sim_df = sim_df[sim_df['Mut_peptide'] != sim_df['WT_peptide']]
    sim_df = sim_df[~sim_df['Mut_peptide'].str.contains('X')]
    
    with open(f"{scratch_dir}/wt_nmer.fasta", "w") as fasta_file:
        for i, (id_val, peptide) in enumerate(zip(sim_df['ID'], sim_df['WT_peptide'])):
            fasta_header = f"{id_val}_{i}"
            fasta_file.write(f">{fasta_header}\n")
            fasta_file.write(f"{peptide}\n")
            
    return f"{scratch_dir}/wt_nmer.fasta"

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


def run_PRIME2(patient_id, processed_alleles, trimmed_alleles,
               similarity_file, mutant_results_file, output_file, scratch_dir):
    """
    Runs PRIME2 once per allele on all WT peptides, merges the WT predictions with
    the mutant predictions and saves the results (with DAI) to a .csv file.

    Parameters:
        patient_id (str): The patient identifier
        processed_alleles (list): A list of trimmed and processed HLA alleles
        trimmed_alleles (list): A list of trimmed HLA alleles
        similarity_file (str): A path to the file mapping WT peptides to mutant peptides
        mutant_results_file (str): A path to the mutant PRIME2 results .csv file
        output_file (str): A path to the output .csv file
        scratch_dir (str): A path to the directory holding the temporary files
    """
    peptides_file = os.path.join(scratch_dir, PEPTIDES_FILENAME)
    prime_output = os.path.join(scratch_dir, PRIME_OUTPUT_FILENAME)

    wt_results = []
    for allele, trimmed_allele in zip(processed_alleles, trimmed_alleles):
        command = ["PRIME", "-i", peptides_file, "-o", prime_output, "-a", allele]
        try:
            subprocess.run(command, check=True, capture_output=True, text=True)
            df_tmp = pd.read_csv(prime_output, comment='#', sep='\t')
        except subprocess.CalledProcessError as e:
            print(f"Failed to run PRIME2 for allele {trimmed_allele}: {e}\n"
                  f"stdout: {e.stdout}\nstderr: {e.stderr}")
            continue

        df_tmp = df_tmp[['Peptide', '%RankBinding_bestAllele']]
        df_tmp.insert(0, 'HLA_allele', trimmed_allele)
        wt_results.append(df_tmp)

    if wt_results:
        wt_df = pd.concat(wt_results, ignore_index=True)
    else:
        wt_df = pd.DataFrame(columns=['HLA_allele', 'Peptide', '%RankBinding_bestAllele'])

    wt_df = wt_df.rename(columns={'Peptide': 'WT_peptide',
                                  '%RankBinding_bestAllele': 'WT_%RankBinding'})
    
    sim_df = pd.read_csv(similarity_file, sep="\t")
    sim_df = sim_df[['WT_peptide', 'Mut_peptide']]
    sim_df = sim_df.rename(columns={'Mut_peptide': 'Peptide'})
    
    mut_df = pd.read_csv(mutant_results_file, sep="\t")
    
    df = pd.merge(wt_df, sim_df, on="WT_peptide", how="left")
    df = pd.merge(mut_df, df, on=['HLA_allele', 'Peptide'], how="left")
    
    df['DAI_%RankBinding'] = df['%RankBinding'] / df['WT_%RankBinding'].replace(0, np.nan)
    
    df = df[['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', '%Rank', 'Score', '%RankBinding', 'DAI_%RankBinding', 'Prediction_tool']]
    df.drop_duplicates(inplace=True)
    df.to_csv(output_file, sep = '\t', index=False)

def cleanup(scratch_dir):
    """
    Cleans up temporary files in the scratch directory.
    """
    for filename in (PEPTIDES_FILENAME, PRIME_OUTPUT_FILENAME, "wt_nmer.fasta"):
        path = os.path.join(scratch_dir, filename)
        if os.path.exists(path):
            os.remove(path)

def main():
    parser = argparse.ArgumentParser(description="Run PRIME2 predictions")
    parser.add_argument("-p", "--patient_id", required=True, help="Patient ID")
    parser.add_argument("--hla_allele_path", required=True, help="Path to the file containing HLA alleles and quality scores")
    parser.add_argument("--hla_quality_threshold", type=int, required=True, help="Quality threshold for filtering alleles")
    parser.add_argument("--similarity_file", required=True, help="Path to the file with self-similarity results.")
    parser.add_argument("--mut_binding_prediction", required=True, help="Path to the file with binding predictions for mutated peptides.")
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
                writer.writerow(['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', '%Rank', 'Score', '%RankBinding', 'DAI_%RankBinding', 'Prediction_tool']) #write an empty header row
            return 

        fasta_path = similarity_to_fasta(args.similarity_file, args.scratch_dir)
        seq_to_identities = process_peptides(fasta_path, args.scratch_dir)
        run_PRIME2(args.patient_id, processed_alleles, trimmed_alleles,
                   args.similarity_file, args.mut_binding_prediction, 
                   args.output_file, args.scratch_dir)
    finally:
        cleanup(args.scratch_dir)
    

if __name__ == "__main__":
    main()