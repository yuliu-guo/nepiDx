import argparse
import pandas as pd
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Generate FASTA from Mut_peptide in a tab-separated input file."
    )

    parser.add_argument(
        "-i", "--input_path",
        required=True,
        help="Path to the input tab-separated file."
    )
    parser.add_argument(
        "-n", "--nmer_output_path",
        required=True,
        help="Output FASTA file path."
    )

    parser.add_argument(
        "--max_hamming_similarity",
        type=float,
        default=1,
        help="Keep rows with hamming_similarity < this value (default: 1)."
    )
    parser.add_argument(
        "--max_blosum62_2_similarity",
        type=float,
        default=1,
        help="Keep rows with blosum62_2_similarity < this value (default: 1)."
    )

    return parser.parse_args()


def validate_columns(df, required_columns):
    missing = [col for col in required_columns if col not in df.columns]
    if missing:
        raise ValueError(
            f"Missing required column(s): {', '.join(missing)}"
        )


def generate_nmer_fasta(
    input_path,
    nmer_output_path,
    max_hamming_similarity=1,
    max_blosum62_2_similarity=1,
):
    df = pd.read_csv(input_path, sep="\t")

    required_columns = [
        "ID",
        "Mut_peptide",
        "WT_peptide",
        "hamming_similarity",
        "blosum62_2_similarity",
    ]
    validate_columns(df, required_columns)

    filtered_df = df[
        (df["hamming_similarity"] < max_hamming_similarity) &
        (df["hamming_similarity"] > 0) &
        (df["blosum62_2_similarity"] < max_blosum62_2_similarity) &
        (df["blosum62_2_similarity"] > 0) &
        (~df["Mut_peptide"].str.contains("X", regex=False, na=False)) &
        (~df["WT_peptide"].str.contains("X", regex=False, na=False))
    ]

    pep_counts = defaultdict(int)

    with open(nmer_output_path, "w") as fasta_file:
        for _, row in filtered_df.iterrows():
            fasta_id = row["ID"]
            sequence = row["Mut_peptide"]

            if pd.notna(sequence) and str(sequence).strip():
                pep_counts[fasta_id] += 1
                fasta_header = f"{fasta_id}_pep{pep_counts[fasta_id]}"

                fasta_file.write(f">{fasta_header}\n")
                fasta_file.write(f"{sequence}\n")

    print(f"FASTA written to: {nmer_output_path}")
    print(f"Number of peptides written: {len(filtered_df)}")


if __name__ == "__main__":
    args = parse_arguments()

    generate_nmer_fasta(
        input_path=args.input_path,
        nmer_output_path=args.nmer_output_path,
        max_hamming_similarity=args.max_hamming_similarity,
        max_blosum62_2_similarity=args.max_blosum62_2_similarity,
    )