import argparse
import pandas as pd


def load_target_specific(filepath):
    """Load target_specific_file and standardize columns."""
    df = pd.read_csv(filepath, sep="\t")

    df = df.rename(columns={
        "Stripped.Sequence": "sequence",
        "n_runs": "n_runs",
        "n_precursors": "n_precursors",
        "Modified.Sequence": "modified_sequence",
        "Protein.Ids": "protein_ids",
        "Protein.Names": "protein_names",
        "Genes": "genes",
        "Precursor.Normalised.max": "precursor_norm_max",
        "Precursor.Normalised.mean": "precursor_norm_mean",
        "Precursor.Normalised.median": "precursor_norm_median",
        "n_target_proteins_matched": "n_proteins_matched",
    })

    df["type"] = "target"
    df["id"] = [f"tar{i}" for i in range(1, len(df) + 1)]

    keep_cols = [
        "id", "sequence", "n_runs", "n_precursors", "modified_sequence",
        "protein_ids", "protein_names", "genes",
        "precursor_norm_max", "precursor_norm_mean", "precursor_norm_median",
        "n_proteins_matched", "type",
    ]
    return df[keep_cols]


def load_mut_specific(filepath):
    """Load mut_specific_file and standardize columns."""
    df = pd.read_csv(filepath, sep="\t")

    df = df.rename(columns={
        "Stripped.Sequence": "sequence",
        "n_runs": "n_runs",
        "n_precursors": "n_precursors",
        "Modified.Sequence": "modified_sequence",
        "Protein.Ids": "protein_ids",
        "Protein.Names": "protein_names",
        "Genes": "genes",
        "Precursor.Normalised.max": "precursor_norm_max",
        "Precursor.Normalised.mean": "precursor_norm_mean",
        "Precursor.Normalised.median": "precursor_norm_median",
        "n_mutant_proteins_matched": "n_proteins_matched",
    })

    df["type"] = "mutant"
    df["id"] = [f"mut{i}" for i in range(1, len(df) + 1)]

    keep_cols = [
        "id", "sequence", "n_runs", "n_precursors", "modified_sequence",
        "protein_ids", "protein_names", "genes",
        "precursor_norm_max", "precursor_norm_mean", "precursor_norm_median",
        "n_proteins_matched", "type",
    ]
    return df[keep_cols]


def merge_dataframes(target_df, mut_df):
    """Concatenate target-specific rows followed by mutant-specific rows."""
    return pd.concat([target_df, mut_df], ignore_index=True)


def main():
    parser = argparse.ArgumentParser(
        description="Gather and index the target-specific and mutant-specific peptides from DIA-NN"
    )
    parser.add_argument(
        "--target_specific_file", required=True,
        help="Path to the target_specific TSV file"
    )
    parser.add_argument(
        "--mut_specific_file", required=True,
        help="Path to the mut_specific TSV file"
    )
    parser.add_argument(
        "--output_file", required=True,
        help="Path to the output TSV file"
    )
    args = parser.parse_args()

    target_df = load_target_specific(args.target_specific_file)
    mut_df = load_mut_specific(args.mut_specific_file)
    merged_df = merge_dataframes(target_df, mut_df)

    merged_df.to_csv(args.output_file, sep="\t", index=False)


if __name__ == "__main__":
    main()