import re
import os
import csv
import pandas as pd
import subprocess
import argparse
from collections import defaultdict, namedtuple

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
    #print(f'full allele is {full_allele}')
    match = re.match(r'^(HLA-[A-Z]\*)(\S+)', full_allele)
    #print(f'match is {match}')
    if not match:
        return full_allele  # fallback — return full allele even if pattern is weird

    prefix, allele_part = match.groups()       # e.g., 'HLA-A*', '01:12'
    pgroup_key = f'{prefix[4:]}{allele_part}'  # strip 'HLA-' -> 'A*01:12'
    
    #print(f"Generated pgroup_map: {pgroup_map}")
    if pgroup_key in pgroup_map:
        #print(f"Generated pgroup_key: {pgroup_key}") 
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

def extract_and_trim_alleles(hla_allele_path, p_group_file, quality_threshold=0):
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
                if quality >= quality_threshold:
                    trimmed = map_to_p_group(allele, pgroup_map)
                    if trimmed:
                        trimmed = trimmed.replace('*', '')  # Remove the asterisk
                        trimmed_alleles.append(trimmed)

    return trimmed_alleles


def run_command_and_log(netmhcpan_cmd, log_file_path):
    # Open the log file in append mode to capture all outputs
    with open(log_file_path, 'a') as log_file:
        try:
            # Run the command and redirect stdout and stderr to the log file
            result = subprocess.run(
                netmhcpan_cmd,
                shell=True,
                check=True,
                stdout=log_file,
                stderr=log_file,
                executable='/bin/bash'
            )
        except subprocess.CalledProcessError as e:
            # Write the error to the log file
            log_file.write(f"Command failed with return code {e.returncode}\n")
            log_file.write(f"Output:\n{e.output}\n")
            log_file.write(f"Error:\n{e.stderr}\n")
            raise

def run_prediction_fasta_and_build_netmhc(patient_id, trimmed_alleles, mutated_peptides_fasta, affinity='YES'):
    """
    Function to load netMHCpan module, run predictions, and parse results directly into a dictionary.

    Parameters:
        patient_id (str): The patient identifier.
        trimmed_alleles (list): A list of trimmed HLA alleles (e.g., HLA-A*02:05).
        mutated_peptides_fasta (str): Path to the mutated peptides FASTA file.
        affinity (str): Whether to include affinity predictions, default is 'YES'.

    Returns:
        net_mhc (dict): A dictionary of parsed netMHCpan results.
    """
    peptide_length = get_first_peptide_length(mutated_peptides_fasta)
    
    hla_pseudo_list_path = "/ref/NetMHCpan4.1_MHC_pseudo.dat"
    
    net_mhc = defaultdict(dict)
    # Mass-Spectrometry Eluted Ligands (EL), Binding Affinity (BA) 
    NetMHCInfo = namedtuple('NetMHCInfo', ['identity', 'score_EL', 'rank_EL', 'score_BA', 'rank_BA', 'affinity'])

    af_col = 15 if affinity == 'YES' else None
    r_col = 16 if affinity == 'YES' else 15


    with open(hla_pseudo_list_path, 'r') as hla_pseudo_list:
        pseudo_alleles = {line.strip().split()[0] for line in hla_pseudo_list}

        for allele in trimmed_alleles:
            #print(f'allele: {allele}')
            #print(type(allele))
            #print("test")
            if allele not in pseudo_alleles:
                #print(pseudo_alleles)
                print(f"WARNING: {allele} cannot be found in hla_pseudo list {hla_pseudo_list_path}")
            else:
                netmhcpan_cmd = (
                    f"module load netMHCpan/4.1b && "
                    f"netMHCpan -a {allele} -f {mutated_peptides_fasta} -l {peptide_length} -BA"
                )
                #TODO: update to netmhcpan/4.1b - this version doesn't recognize * in alleles ( HLA-A*02:01 should be cleaned as HLA-A02:01)

                try:
                    # Run the command and capture the output
                    result = subprocess.run(netmhcpan_cmd, shell=True, check=True, capture_output=True, text=True)
                    output = result.stdout
                    

                    # Parse the output directly into the net_mhc dictionary
                    for line in output.splitlines():
                        if re.search(r'^\s+\d+\s', line):
                            print(f'line is {line}')
                            line = [x.strip() for x in line.split()]
                            line[2] = '-' * len(line[2]) if line[2].startswith('XXXX') else line[2]
                            HLA_allele = line[1].replace('*', '')
                            peptide = line[2]
                            identity = line[10]
                            score_EL = float(line[11])
                            rank_EL = float(line[12])
                            score_BA = float(line[13])
                            rank_BA = float(line[14])
                            affinity = float(line[af_col]) if af_col is not None else None
                            netmhc_info = NetMHCInfo(identity, score_EL, rank_EL, score_BA, rank_BA, affinity)
                            net_mhc[HLA_allele][peptide] = netmhc_info


                    print(f"Successfully ran and parsed netMHCpan for allele {allele}")
                    
                except subprocess.CalledProcessError as e:
                    print(f"Failed to run netMHCpan for allele {allele}: {e}")
                    continue

    #print(f"Completed netMHCpan prediction for patient: {patient_id}. ")
    return net_mhc


def write_net_mhc_to_csv(patient_id, net_mhc, output_file):
    """
    Write the netMHC dictionary data to an csv file.

    Parameters:
        patient_id (str): The patient identifier.
        net_mhc (dict): The netMHC prediction results.
        output_file (str): Directory to save the csv file.
    """
    # if net_mhc is None:
    #     # Create an empty file
    #     open(output_file, 'w').close()
    #     print(f"netMHC data is None for patient: {patient_id}. Empty file saved to {output_file}")
    # else:
    data = []
    prediction_tool = 'netMHCpan'
    for HLA_allele, peptides in net_mhc.items():
        for peptide, info in peptides.items():
            data.append([
                patient_id, 
                info.identity,
                HLA_allele,
                peptide,
                info.score_EL,
                info.rank_EL,
                info.score_BA,
                info.rank_BA,
                info.affinity,
                prediction_tool
            ])
    
    # Include 'Patient_ID' as a new column
    df = pd.DataFrame(data, columns=['Patient_ID', 'Identity', 'HLA_allele', 'Peptide', 'Score_EL', 'Rank_EL', 'Score_BA', 'Rank_BA', 'Affinity', 'Prediction_tool'])

    df.to_csv(output_file, sep = '\t', index=False)
    print(f"Completed netMHCpan prediction for patient: {patient_id}. Output saved to {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Run netMHCpan predictions')
    parser.add_argument("-p", "--patient_id", required=True, help="Patient ID")
    parser.add_argument("--hla_allele_path", required=True, help="Path to the file containing HLA alleles and quality scores")
    parser.add_argument("--hla_quality_threshold", type=int, required=True, help="Quality threshold for filtering alleles")
    parser.add_argument("--mutated_peptides_fasta", required=True, help="Path to the mutated peptides FASTA file")
    parser.add_argument("-o", "--output_file", required=True, help="Path to save the output files")
    parser.add_argument("--p_group_file", type=str, required=True, help="Path to P_GROUP_FILE for HLA allele mapping")

    # Parse arguments
    args = parser.parse_args()

    # Extract and trim HLA alleles
    trimmed_alleles = extract_and_trim_alleles(args.hla_allele_path, args.p_group_file, args.hla_quality_threshold)
    if not trimmed_alleles:
        print("No allele passed the quality threshold")
        with open(args.output_file, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow([
                  'Patient_ID', 'Identity', 'HLA_allele', 'Peptide', 'Score_EL', 
                  'Rank_EL', 'Score_BA', 'Rank_BA', 'Affinity', 'Prediction_tool'
            ]) #write an empty header row
        return 

    # Run netMHCpan predictions and directly parse the results into a dictionary
    #net_mhc = run_prediction_fasta_and_build_netmhc(args.patient_id, trimmed_alleles, args.mutated_peptides_fasta, args.peptide_length)

    # Write the results to an csv file
    write_net_mhc_to_csv(
        args.patient_id, 
        run_prediction_fasta_and_build_netmhc(
            args.patient_id, trimmed_alleles, args.mutated_peptides_fasta
            ), 
            args.output_file
        )


if __name__ == "__main__":
    main()


