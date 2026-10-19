import os
import ast
import pandas as pd
import argparse

def parse_arguments():
    """
    Parse command-line arguments for generating FASTA files.
    """
    parser = argparse.ArgumentParser(description="Generate FASTA files for mutated peptides.")

    parser.add_argument(
        "-i", "--input_path",
        type=str,
        required=True,
        help="Path to the input file."
    )
    parser.add_argument(
        "-n", "--nmer_output_path",
        type=str,
        required=False,
        help="Output file path for 'nmer' FASTA."
    )
    parser.add_argument(
        "-f", "--flanked_output_path",
        type=str,
        required=False,
        help="Output file path for 'flanked_mutated_protein_sequence' FASTA."
    )
    parser.add_argument(
        "-m", "--min_alt_read_base",
        type=int,
        default=0,
        help="Minimum alt_read_base_counts for filtering rows (default: 0)."
    )
    parser.add_argument(
        "-e", "--expand_peptide_output_path",
        type=str,
        required=False,
        help="Output file path for dataframe that is expanded with peptide."
    )
    

    return parser.parse_args()



def generate_nmer_fasta(df, nmer_output_path, min_alt_read_base=0):
    """
    Generate a FASTA file for peptides from the 'nmer' column.
    
    Arguments:
    df (DataFrame): Input DataFrame containing mutation data.
    nmer_output_path (str): File path to save FASTA for 'nmer' peptides.
    min_alt_read_base (int): Minimum RNA_ALT_N for filtering rows.
    """
    filtered_df = df[df['RNA_ALT_N'] >= min_alt_read_base]
    
    with open(nmer_output_path, 'w') as nmer_fasta_file:
        for index, row in filtered_df.iterrows():
            nmer_str = row.get('nmer', None)
            fasta_id = row.get('ID', '')
            #chrom = row.get('CHROM', '')
            #pos = row.get('POS', '')
            #symbol = row.get('SYMBOL', '')
            
            if pd.notna(nmer_str) and nmer_str:
                try:
                    # Evaluate the string into a list of peptides
                    nmer_list = ast.literal_eval(nmer_str)
                    for peptide in nmer_list:
                        #fasta_id = f"{chrom}:{pos}:{symbol}"
                        #fasta_id = f"{chrom}_{pos}"
                        
                        nmer_fasta_file.write(f">{fasta_id}\n")
                        nmer_fasta_file.write(f"{peptide}\n")
                except Exception as e:
                    print(f"Error processing 'nmer' in row {index}: {e}")

    print(f"FASTA file for 'nmer' peptides written to: {nmer_output_path}")


def generate_flanked_fasta(df, flanked_output_path, min_alt_read_base=0):
    """
    Generate a FASTA file for sequences from the 'flanked_mutated_protein_sequence' column.
    
    Arguments:
    df (DataFrame): Input DataFrame containing mutation data.
    flanked_output_path (str): File path to save FASTA for 'flanked_mutated_protein_sequence'.
    min_alt_read_base (int): Minimum RNA_ALT_N for filtering rows.
    """

    filtered_df = df[df['RNA_ALT_N'] >= min_alt_read_base]

    with open(flanked_output_path, 'w') as flanked_fasta_file:
        for index, row in filtered_df.iterrows():
            #fasta_id = row.get('ID', '')
            #chrom = row.get('CHROM', '')
            #pos = row.get('POS', '')
            #symbol = row.get('SYMBOL', '')
            sequence = row.get('flanked_mutated_protein_sequence', None)
            
            if pd.notna(sequence) and sequence:
                #fasta_id = f"{chrom}:{pos}:{symbol}"
                #fasta_id = f"{chrom}_{pos}"
                fasta_id = row.get('ID', '')
                flanked_fasta_file.write(f">{fasta_id}\n")
                flanked_fasta_file.write(f"{sequence}\n")
    
    print(f"FASTA file for 'flanked_mutated_protein_sequence' written to: {flanked_output_path}")


def generate_fasta_files(input_path, nmer_output_path=None, flanked_output_path=None, expand_peptide_output_path=None, min_alt_read_base=0):
    """
    Generate FASTA files based on the specified output paths.
    
    Arguments:
    input_path (str): Path to the input file containing mutation data.
    nmer_output_path (str): File path to save FASTA for 'nmer' peptides (optional).
    flanked_output_path (str): File path to save FASTA for 'flanked_mutated_protein_sequence' (optional).
    min_alt_read_base (int): Minimum alt_read_base_counts for filtering rows.
    """
    df = pd.read_csv(input_path,  sep='\t')
    
    if nmer_output_path:
        generate_nmer_fasta(df, nmer_output_path, min_alt_read_base)
    
    if flanked_output_path:
        generate_flanked_fasta(df, flanked_output_path, min_alt_read_base)

    if expand_peptide_output_path:
        expand_df_nmer(df, expand_peptide_output_path, min_alt_read_base)

    if not nmer_output_path and not flanked_output_path and not expand_peptide_output_path:
        print("No output paths provided. Nothing to generate.")


import ast
def expand_peptide_df(df, expand_peptide_output_path, min_alt_read_base=0):
    """
    expand summary with each nmer peptide is one row
    Combining with Chrom, Pos, SYMBOL, alt_read

    """
    print("Columns in the input file:", df.columns)

    columns = ['CHROM','POS','SYMBOL','HGVSp','gene_count_TPM','alt_read_base_counts','nmer']
    df1 = pd.DataFrame(df, columns=columns)
    df1.rename(columns={'nmer':'Peptide'}, inplace = True)
    # Convert the 'nmer' column to lists if necessary
    df1['Peptide'] = df1['Peptide'].apply(lambda x: ast.literal_eval(x) if isinstance(x, str) else x)
    df_filtered = df1[df1['alt_read_base_counts'] >= min_alt_read_base]
    # Explode the 'nmer' column
    df_expanded = df_filtered.explode('Peptide').reset_index(drop=True)

    # Write DataFrame to CSV
    df_expanded.to_csv(expand_peptide_output_path, sep=';', index=False)

    return df_expanded


if __name__ == "__main__":
    # Parse arguments from the command line
    args = parse_arguments()
    
#    Call the function to generate FASTA files
    generate_fasta_files(
        input_path=args.input_path,
        nmer_output_path=args.nmer_output_path,
        flanked_output_path=args.flanked_output_path,
        expand_peptide_output_path=args.expand_peptide_output_path,
        min_alt_read_base=args.min_alt_read_base
    )

