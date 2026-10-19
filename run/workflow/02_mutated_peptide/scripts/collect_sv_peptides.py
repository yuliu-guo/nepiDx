#!/usr/bin/env python3
import os
import argparse
import pandas as pd


def parse_arguments():
    """
    Parse command-line arguments for generating FASTA files.
    """
    parser = argparse.ArgumentParser(
        description="Generate a FASTA file from TSV column 'novel_kmers_sequence' (semicolon-separated)."
    )
    parser.add_argument(
        "-i", "--input_path",
        type=str,
        required=True,
        help="Path to the input TSV file (must contain columns: ID, novel_kmers_sequence)."
    )
    parser.add_argument(
        "-o", "--nmer_output_path",
        type=str,
        required=False,
        default=None,
        help="Output FASTA path. Default: <input_basename>.novel_kmers.fasta"
    )
    return parser.parse_args()


def extract_peptides_semicolon(cell) -> list[str]:
    """
    novel_kmers_sequence is a string like: 'AAA;BBB;CCC'
    Returns ['AAA','BBB','CCC'].
    """
    if cell is None:
        return []
    s = str(cell).strip()
    if not s or s.lower() == "nan":
        return []
    return [p.strip() for p in s.split(";") if p.strip()]


def generate_nmer_fasta(df: pd.DataFrame, nmer_output_path: str):
    """
    Writes FASTA records as:
      >ID_pep1
      PEPTIDE
      >ID_pep2
      PEPTIDE
      ...
    """
    required = {"ID", "novel_kmers_sequence"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Input file is missing required column(s): {', '.join(sorted(missing))}")

    n_written = 0
    n_rows_with_any = 0

    with open(nmer_output_path, "w") as f:
        for idx, row in df.iterrows():
            fasta_id = str(row["ID"]).strip()  # you said it will never be empty

            peptides = extract_peptides_semicolon(row.get("novel_kmers_sequence", None))
            if peptides:
                n_rows_with_any += 1

            for j, pep in enumerate(peptides, start=1):
                f.write(f">{fasta_id}_pep{j}\n{pep}\n")
                n_written += 1

    print(f"Wrote {n_written} FASTA records from {n_rows_with_any} rows to: {nmer_output_path}")


def main():
    args = parse_arguments()

    if not os.path.exists(args.input_path):
        raise FileNotFoundError(f"Input file not found: {args.input_path}")

    out_path = args.nmer_output_path
    if out_path is None:
        base, _ = os.path.splitext(args.input_path)
        out_path = base + ".novel_kmers.fasta"

    # Your pipeline writes TSV
    df = pd.read_csv(args.input_path, sep="\t", dtype=str)

    generate_nmer_fasta(df, out_path)


if __name__ == "__main__":
    main()

