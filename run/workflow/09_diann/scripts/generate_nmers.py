import argparse
import pandas as pd


def load_input_file(filepath):
    """Load the merged peptide TSV file."""
    return pd.read_csv(filepath, sep="\t")


def generate_nmers(df, nmer_length):
    """Generate all nmers of a given length from each sequence.

    Rows whose sequence is shorter than nmer_length are skipped entirely.
    Returns a list of (header, nmer_sequence) tuples.
    """
    records = []

    for _, row in df.iterrows():
        seq_id = row["id"]
        sequence = row["sequence"]

        if len(sequence) < nmer_length:
            continue

        n_nmers = len(sequence) - nmer_length + 1
        for i in range(1, n_nmers + 1):
            nmer = sequence[i - 1:i - 1 + nmer_length]
            header = f">{seq_id}_pep{i}"
            records.append((header, nmer))

    return records


def write_fasta(records, output_file):
    """Write (header, sequence) records to a FASTA file."""
    with open(output_file, "w") as f:
        for header, sequence in records:
            f.write(f"{header}\n{sequence}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Generate nmer FASTA file from target/mutant-specific peptides"
    )
    parser.add_argument(
        "--input_file", required=True,
        help="Path to the input TSV file (output of the merge script)"
    )
    parser.add_argument(
        "--nmer_length", required=True, type=int,
        help="Length of the nmers to generate"
    )
    parser.add_argument(
        "--output_file", required=True,
        help="Path to the output FASTA file"
    )
    args = parser.parse_args()

    df = load_input_file(args.input_file)
    records = generate_nmers(df, args.nmer_length)
    write_fasta(records, args.output_file)


if __name__ == "__main__":
    main()