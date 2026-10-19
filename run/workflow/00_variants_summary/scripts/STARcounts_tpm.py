#!/usr/bin/env python3

import os
import warnings
import numpy as np
import pandas as pd

# gene length table path
gene_length_table_path = "/ref/gene_length_table.txt"


# Internal caches 
_STAR_TPM_LOOKUP_CACHE = {}
_GENE_LENGTH_DF_CACHE = {}


# File readers / processors
def read_STARcounts(STARcounts_path):
    """
    Read the STARcounts file and return a dataframe with columns: ['gene', 'STDIN'].

    Accepts either:
      1) 2-column format: gene \\t counts
      2) 4-column format: gene \\t unstrand \\t strand1 \\t strand2 (STAR --quantMode GeneCounts)
         -> maps 'unstrand' to 'STDIN' to preserve downstream logic.

    If the path is missing/invalid, returns an empty DataFrame with the expected columns.
    """
    if not STARcounts_path or str(STARcounts_path).strip() in ("", ".") or not os.path.exists(STARcounts_path):
        warnings.warn(
            f"STARcounts_path is missing or invalid: {STARcounts_path!r}. Returning empty DataFrame.",
            UserWarning
        )
        return pd.DataFrame(columns=["gene", "STDIN"])

    df_raw = pd.read_csv(
        STARcounts_path,
        sep="\t",
        header=None,
        comment="#",
        dtype=str,
        engine="python"
    )

    if df_raw.empty:
        return pd.DataFrame(columns=["gene", "STDIN"])

    # Use only the first 4 columns if extra columns exist
    if df_raw.shape[1] > 4:
        df_raw = df_raw.iloc[:, :4]

    def _drop_header_if_present(frame, expected_first_col_names):
        first = frame.iloc[0].astype(str).str.lower().tolist()
        if first and (
            first[0] in expected_first_col_names
            or any(x in ("stdin", "count", "counts", "unstrand", "unstranded", "strand1", "strand2") for x in first[1:])
        ):
            return frame.iloc[1:].reset_index(drop=True)
        return frame

    ncols = df_raw.shape[1]

    if ncols == 2:
        df = _drop_header_if_present(df_raw, expected_first_col_names={"gene", "geneid", "symbol"})
        df.columns = ["gene", "STDIN"]
        return df

    if ncols >= 4:
        df4 = _drop_header_if_present(df_raw.iloc[:, :4], expected_first_col_names={"gene", "geneid", "symbol"})
        df4.columns = ["gene", "unstrand", "strand1", "strand2"]
        df = df4[["gene", "unstrand"]].copy()
        df.columns = ["gene", "STDIN"]
        return df

    raise ValueError(f"Unrecognized STAR counts layout with {ncols} columns in file: {STARcounts_path}")


def filter_N_rows(df):
    """
    Filter out rows with N_, such as N_unmapped, N_multimapping, N_noFeature, N_ambiguous
    """
    if df.empty:
        return df.copy()
    return df[~df["gene"].astype(str).str.startswith("N_")].copy()


def read_gene_length_table(path):
    """
    Read the gene_length_table.txt and return a dataframe with selected columns.
    Cached so it's read once.
    """
    if path in _GENE_LENGTH_DF_CACHE:
        return _GENE_LENGTH_DF_CACHE[path]

    gene_length_df = pd.read_csv(path, sep="\t")
    gene_length_df = gene_length_df[["gene_name", "gene_length"]].copy()

    _GENE_LENGTH_DF_CACHE[path] = gene_length_df
    return gene_length_df


def merge_STARcounts_with_gene_length(STARcounts_df, gene_length_df):
    """
    Merge the STARcounts DataFrame with gene length table on gene name
    """
    if STARcounts_df.empty:
        out = STARcounts_df.copy()
        out["gene_length"] = np.nan
        return out

    merged_df = STARcounts_df.merge(gene_length_df, left_on="gene", right_on="gene_name", how="left")
    merged_df = merged_df.drop(columns=["gene_name"])
    return merged_df


def calculate_tpm(df, counts_column="STDIN"):
    """
    Calculates TPM for a specific sample column.
    Step 1: Normalize for gene length (kilobases) -> Reads Per Kilobase (RPK)
    Step 2: Normalize for sequencing depth (divide by total RPK)
    """
    if df.empty:
        df = df.copy()
        df["TPM"] = []
        return df

    df = df.copy()
    df[counts_column] = pd.to_numeric(df[counts_column], errors="coerce").fillna(0.0)
    df["gene_length"] = pd.to_numeric(df["gene_length"], errors="coerce")

    # Convert gene length to kilobases
    df["gene_length_kb"] = df["gene_length"] / 1000.0

    # RPK
    df["RPK"] = df[counts_column] / df["gene_length_kb"]

    # Scaling factor
    scaling_factor = df["RPK"].replace([np.inf, -np.inf], np.nan).sum(skipna=True)

    if not np.isfinite(scaling_factor) or scaling_factor == 0:
        warnings.warn(
            "TPM scaling factor is 0 or non-finite (likely missing/zero gene lengths). Setting all TPM to 0.0.",
            UserWarning
        )
        df["TPM"] = 0.0
        return df

    df["TPM"] = (df["RPK"] / scaling_factor) * 1e6
    df["TPM"] = df["TPM"].replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return df


def get_gene_tpm(df, gene_name):
    """
    Return TPM value for the specified gene (DataFrame path; kept for compatibility).
    """
    if df.empty:
        warnings.warn(f"Warning: gene '{gene_name}' not found in empty dataset. Returning 0.0.", UserWarning)
        return 0.0

    row = df[df["gene"] == gene_name]
    if row.empty:
        warnings.warn(f"Warning: gene '{gene_name}' not found in dataset. Returning 0.0.", UserWarning)
        return 0.0

    tpm_value = float(row["TPM"].values[0])
    if not np.isfinite(tpm_value):
        warnings.warn(f"Warning: TPM for '{gene_name}' is non-finite (nan/inf). Returning 0.0.", UserWarning)
        return 0.0
    return tpm_value



# load once -> dict lookup
def load_star_tpm_lookup(STARcounts_path):
    """
    Load STARcounts once, compute TPM for ALL genes once, return lookup dict: gene -> TPM.
    Cached by STARcounts_path.
    """
    if STARcounts_path in _STAR_TPM_LOOKUP_CACHE:
        return _STAR_TPM_LOOKUP_CACHE[STARcounts_path]

    df = read_STARcounts(STARcounts_path)
    df_filtered = filter_N_rows(df)

    gene_length_df = read_gene_length_table(gene_length_table_path)

    df_merged = merge_STARcounts_with_gene_length(df_filtered, gene_length_df)
    df_normalized = calculate_tpm(df_merged)

    # Build lookup dict; use 0.0 for non-finite values
    tpm_series = pd.to_numeric(df_normalized["TPM"], errors="coerce").replace([np.inf, -np.inf], np.nan).fillna(0.0)
    lookup = dict(zip(df_normalized["gene"].astype(str), tpm_series.astype(float)))

    _STAR_TPM_LOOKUP_CACHE[STARcounts_path] = lookup
    return lookup



# Backward compatible main()
def main(STARcounts_path, gene_name):
    """
    Backward-compatible function:
    - Now uses the cached lookup (fast for repeated calls).
    - Returns 0.0 if gene not present (same behavior as get_gene_tpm warnings path).
    """
    lookup = load_star_tpm_lookup(STARcounts_path)
    val = lookup.get(gene_name)
    if val is None:
        warnings.warn(f"Warning: gene '{gene_name}' not found in dataset. Returning 0.0.", UserWarning)
        return 0.0
    if not np.isfinite(val):
        warnings.warn(f"Warning: TPM for '{gene_name}' is non-finite (nan/inf). Returning 0.0.", UserWarning)
        return 0.0
    return float(val)