import os
import subprocess
import argparse
import pandas as pd
import numpy as np
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
                        processed = trimmed.replace('*', '')  # Remove the asterisk
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

def _write_empty_netmhcstabpan_tsv(out_path: str, include_affinity: bool) -> None:
    """
    Create an empty netMHCstabpan TSV that is compatible with:
      pd.read_csv(out_path, skiprows=1, sep='\t')
    i.e. first line is a comment, second line is the header.
    """
    # These are the columns netMHCstabpan typically writes in the -xls TSV,
    # including 'Pos' which your code drops after reading.
    netmhc_cols = ['Pos', 'Peptide', 'Identity', 'Pred', 'Thalf(h)', '%Rank_stab']
    if include_affinity:
        netmhc_cols.extend(['1-log50k', 'Aff(nM)', '%Rank_aff', 'Comb', '%Rank_comb'])
    netmhc_cols.extend(['Ave', 'NB'])

    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("# netMHCstabpan output (empty)\n")
        f.write("\t".join(netmhc_cols) + "\n")


def run_netMHCstabpan(patient_id, processed_alleles, trimmed_alleles, similarity_file, mutant_results_file, output_file, scratch_dir, include_affinity=True):
    """
    Runs netMHCstabpan prediction with given alleles and peptides. Saves the result to a .csv file

    Parameters:
        patient_id (str): The patient identifier
        processed_alleles (list): A list of trimmed and processed HLA alleles
        trimmed_alleles (list): A list of trimmed HLA alleles
        peptides_file (str): A path to fasta file with mutated peptides
        output_file (str): A path to the output .csv file
        scratch_dir (str): A directory for storing temporary files; the files are removed after the process finishes
        include_affinity (boolean): Whether to include affinity predictions ("-ia" flag in netMHCstabpan); default: True
    """
    columns = ['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Pred', 'Thalf(h)', '%Rank_stab']
    if include_affinity:
        columns.extend(['1-log50k', 'Aff(nM)', '%Rank_aff', 'Comb', '%Rank_comb'])
    columns.extend(['Ave', 'NB', 'Prediction_tool'])
    wt_df = pd.DataFrame(columns=columns)
    
    os.makedirs(scratch_dir, exist_ok=True)
    out_path = os.path.join(scratch_dir, "netMHCstabpan_output.tsv")
    
    peptides_file = similarity_to_fasta(similarity_file, scratch_dir)
    peptide_length = get_first_peptide_length(peptides_file)
    
    # For some samples, there is no peptide that passed filter on min_alt_count = 1, therefore {peptides_file} is empty(or could be there is no valid alleles) -> 
    # write empty outputs and return will save the pipeline from crash. 
    peptides_missing_or_empty = (
        (not peptides_file) or
        (not os.path.exists(peptides_file)) or
        (os.path.getsize(peptides_file) == 0)
    )
    if peptides_missing_or_empty or (not processed_alleles):
        _write_empty_netmhcstabpan_tsv(out_path, include_affinity=include_affinity)
        wt_df.to_csv(output_file, sep='\t', index=False)
        return wt_df


    for i, allele in enumerate(processed_alleles):
        # Skip empty/None allele safely (no subprocess call)
        if (allele is None) or (not str(allele).strip()):
            continue
            
        if include_affinity:
            command = (
                f"netMHCstabpan -ia -f {peptides_file} -a {allele} -l {peptide_length} -xls -xlsfile {out_path}"
            )
        else:
            command = (
                f"netMHCstabpan -f {peptides_file} -a {allele} -l {peptide_length} -xls -xlsfile {out_path}"
            )
        try:
            result = subprocess.run(command, capture_output=True, check=True, text=True, shell=True)
            
            # If netMHCstabpan didn't produce the file, create an empty one so parsing won't crash
            if not os.path.exists(out_path):
                _write_empty_netmhcstabpan_tsv(out_path, include_affinity=include_affinity)

                print("Expected output missing:", out_path)
                print("scratch_dir files:", os.listdir(scratch_dir))
                print("STDERR (first 2000):", result.stderr[:2000])
                print("STDOUT (first 2000):", result.stdout[:2000])

            #####################
            df_tmp = pd.read_csv(out_path, skiprows=1, sep='\t')
            df_tmp = df_tmp.drop(columns=['Pos'])
            df_tmp.insert(0, 'Patient_ID', patient_id)
            df_tmp.insert(1, 'HLA_allele', trimmed_alleles[i])
            df_tmp['Prediction_tool'] = 'netMHCstabpan'
            df_tmp.columns = wt_df.columns
            wt_df = pd.concat([wt_df, df_tmp], ignore_index=True)
        except subprocess.CalledProcessError as e:
            print(f"Failed to run netMHCstabpan for allele {trimmed_alleles[i]}: {e}")
            continue
    
    wt_columns = ['HLA_allele', 'Peptide', '%Rank_stab']
    if include_affinity:
        wt_columns.extend(['%Rank_comb'])
        
    wt_df = wt_df[wt_columns]
    wt_df = wt_df.rename(columns={'Peptide': 'WT_peptide',
                                  '%Rank_stab': 'WT_%Rank_stab'})
    if include_affinity:
        wt_df = wt_df.rename(columns={'%Rank_comb': 'WT_%Rank_comb'})
                                  
    sim_df = pd.read_csv(similarity_file, sep="\t")
    sim_df = sim_df[['WT_peptide', 'Mut_peptide']]
    sim_df = sim_df.rename(columns={'Mut_peptide': 'Peptide'})
    
    mut_df = pd.read_csv(mutant_results_file, sep="\t")
    
    df = pd.merge(wt_df, sim_df, on="WT_peptide", how="left")
    df = pd.merge(mut_df, df, on=['HLA_allele', 'Peptide'], how="left")
    
    df['DAI_%Rank_stab'] = df['%Rank_stab'] / df['WT_%Rank_stab'].replace(0, np.nan)
    if include_affinity:
        df['DAI_%Rank_comb'] = df['%Rank_comb'] / df['WT_%Rank_comb'].replace(0, np.nan)
    
    result_columns = ['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Pred', 'Thalf(h)', '%Rank_stab']
    if include_affinity:
        result_columns.extend(['1-log50k', 'Aff(nM)', '%Rank_aff', 'Comb', '%Rank_comb'])
    result_columns.extend(['DAI_%Rank_stab'])
    if include_affinity:
        result_columns.extend(['DAI_%Rank_comb'])
    result_columns.extend(['Prediction_tool'])
        
    df = df[result_columns]
    df.drop_duplicates(inplace=True)
    df.to_csv(output_file, sep = '\t', index=False)    

    

def cleanup(scratch_dir):
    """
    Cleans up temporary files in the scratch directory.
    """
    if os.path.exists(os.path.join(scratch_dir, "netMHCstabpan_output.tsv")):
        os.remove(os.path.join(scratch_dir, "netMHCstabpan_output.tsv"))

def main():
    parser = argparse.ArgumentParser(description="Run netMHCstabpan predictions")
    parser.add_argument("-p", "--patient_id", required=True, help="Patient ID")
    parser.add_argument("--hla_allele_path", required=True, help="Path to the file containing HLA alleles and quality scores")
    parser.add_argument("--hla_quality_threshold", type=int, required=True, help="Quality threshold for filtering alleles")
    parser.add_argument("--similarity_file", required=True, help="Path to the file with self-similarity results.")
    parser.add_argument("--mut_binding_prediction", required=True, help="Path to the file with binding predictions for mutated peptides.")
    parser.add_argument("-o", "--output_file", required=True, help="Path to save the output files")
    parser.add_argument("--scratch_dir", type=str, required=True, help="Path to the scratch directory")
    parser.add_argument("-ia", "--include_affinity", action="store_true", help="Specifies whether to include affinity predictions")
    parser.add_argument("--p_group_file", type=str, required=True, help="Path to P_GROUP_FILE for HLA allele mapping")

    args = parser.parse_args()

    try:
        processed_alleles, trimmed_alleles = extract_and_trim_alleles(args.hla_allele_path, args.p_group_file, args.hla_quality_threshold)
        if not trimmed_alleles:
            print("No allele passed the quality threshold")
            columns = ['Patient_ID', 'HLA_allele', 'Peptide', 'Identity', 'Pred', 'Thalf(h)', '%Rank_stab']
            if args.include_affinity:
                columns.extend(['1-log50k', 'Aff(nM)', '%Rank_aff', 'Comb', '%Rank_comb'])
            columns.extend(['DAI_%Rank_stab'])
            if args.include_affinity:
                columns.extend(['DAI_%Rank_comb'])
            columns.extend(['Prediction_tool'])
            
            with open(args.output_file, 'w', newline='') as csvfile:
                writer = csv.writer(csvfile, delimiter='\t')
                writer.writerow(columns) #write an empty header row
            return 

        run_netMHCstabpan(args.patient_id, processed_alleles, trimmed_alleles,
                            args.similarity_file, args.mut_binding_prediction, args.output_file,
                            args.scratch_dir, args.include_affinity)
    finally:
        cleanup(args.scratch_dir)
    

if __name__ == "__main__":
    main()