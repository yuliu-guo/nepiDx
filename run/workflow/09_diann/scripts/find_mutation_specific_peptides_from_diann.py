#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Find mutation-specific peptides from a DIA-NN output folder.

Usage:
    python find_mutation_specific_peptides_from_diann.py /path/to/diann_folder
    python find_mutation_specific_peptides_from_diann.py /path/to/diann_folder --mutant-prefix MUTANT
    python find_mutation_specific_peptides_from_diann.py /path/to/diann_folder --outdir custom_results

Requirements:
    pip install pandas pyarrow
"""

from pathlib import Path
import pandas as pd
import re
import argparse
import sys


DEFAULT_MUTANT_PREFIX = "MUTANT"
DEFAULT_MAX_QVALUE = 0.01
DEFAULT_MAX_GLOBAL_QVALUE = 0.01

REPORT_SEARCH_COLUMNS = [
    "Protein.Ids",
    "Protein.Group",
    "Protein.Names",
    "Genes",
    "Gene.Group",
]

REPORT_KEEP_COLUMNS = [
    "Run",
    "Precursor.Id",
    "Modified.Sequence",
    "Stripped.Sequence",
    "Precursor.Charge",
    "Protein.Ids",
    "Protein.Group",
    "Protein.Names",
    "Genes",
    "Gene.Group",
    "Q.Value",
    "Global.Q.Value",
    "Lib.Q.Value",
    "PEP",
    "Precursor.Quantity",
    "Precursor.Normalised",
    "RT",
    "Predicted.RT",
    "IM",
    "Predicted.IM",
]

QUANTITY_COLUMNS = [
    "Precursor.Normalised",
    "Precursor.Quantity",
    "Quantity",
]

PROTSEQ_ID_COLUMNS = [
    "Protein.Ids",
    "Protein",
    "Protein.ID",
    "Protein.Id",
    "ProteinID",
    "ProteinIds",
]

PROTSEQ_NAME_COLUMNS = [
    "Protein.Names",
    "Protein.Name",
    "ProteinName",
    "ProteinNames",
]

PROTSEQ_GENE_COLUMNS = [
    "Genes",
    "Gene.Names",
    "Gene.Name",
    "GeneNames",
    "Gene",
]

PROTSEQ_SEQUENCE_COLUMNS = [
    "Sequence",
    "Protein.Sequence",
    "ProteinSequence",
    "Stripped.Sequence",
]



def find_first_existing(df: pd.DataFrame, candidates: list[str]) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    return None


def existing_columns(df: pd.DataFrame, cols: list[str]) -> list[str]:
    return [c for c in cols if c in df.columns]


def choose_quantity_column(df: pd.DataFrame) -> str | None:
    for c in QUANTITY_COLUMNS:
        if c in df.columns:
            return c
    return None


def clean_sequence(seq: str) -> str:
    if pd.isna(seq):
        return ""
    seq = str(seq).upper()
    seq = re.sub(r"[^A-Z]", "", seq)
    return seq


def safe_contains_prefix(value: str, prefix: str) -> bool:
    if pd.isna(value):
        return False
    return str(value).upper().startswith(prefix.upper())


def row_has_mutant_prefix(row: pd.Series, columns: list[str], prefix: str) -> bool:
    for col in columns:
        if col in row.index and safe_contains_prefix(row[col], prefix):
            return True
    return False


def filter_qvalues(df: pd.DataFrame, max_qvalue: float, max_global_qvalue: float) -> pd.DataFrame:
    out = df.copy()
    if "Q.Value" in out.columns:
        out = out[out["Q.Value"].fillna(1.0) <= max_qvalue]
    if "Global.Q.Value" in out.columns:
        out = out[out["Global.Q.Value"].fillna(1.0) <= max_global_qvalue]
    return out


def build_mutant_tag_mask(df: pd.DataFrame, search_cols: list[str], prefix: str) -> pd.Series:
    cols = existing_columns(df, search_cols)
    if not cols:
        raise ValueError(f"None of these columns exist in report.parquet: {search_cols}")

    mask = pd.Series(False, index=df.index)
    for col in cols:
        col_mask = (
            df[col]
            .fillna("")
            .astype(str)
            .str.upper()
            .str.contains(re.escape(prefix.upper()), regex=True)
        )
        mask = mask | col_mask
    return mask


def find_matching_proteins(peptide: str, protein_df: pd.DataFrame, seq_col: str, id_col: str) -> list[str]:
    pep = clean_sequence(peptide)
    if not pep:
        return []

    hits = []
    for _, row in protein_df.iterrows():
        seq = row[seq_col]
        if pep in seq:
            hits.append(str(row[id_col]))
    return hits



def discover_report_file(diann_dir: Path) -> Path:
    """
    Discover the DIA-NN main report parquet.

    Accept:
      - <prefix>.parquet

    Exclude:
      - <prefix>_outlib.parquet
      - <prefix>.site_report.parquet
      - any *.xic.parquet
      - any parquet inside xic/ subfolders
    """
    candidates = []

    for p in diann_dir.glob("*.parquet"):
        name = p.name

        # exclude non-main parquet outputs
        if name.endswith("_outlib.parquet"):
            continue
        if name.endswith(".site_report.parquet"):
            continue
        if name.endswith(".xic.parquet"):
            continue

        candidates.append(p)

    if len(candidates) == 0:
        raise FileNotFoundError(
            f"No DIA-NN main report parquet found in: {diann_dir}"
        )

    if len(candidates) == 1:
        return candidates[0]

    # Prefer the simplest "<prefix>.parquet" file name
    simple_candidates = [
        p for p in candidates
        if p.stem and not p.stem.endswith("_outlib") and not p.stem.endswith(".site_report")
    ]

    if len(simple_candidates) == 1:
        return simple_candidates[0]

    raise FileNotFoundError(
        f"Could not uniquely determine main report parquet in: {diann_dir}\n"
        f"Candidates found: {[p.name for p in candidates]}"
    )

def discover_protein_description_file(diann_dir: Path) -> Path:
    candidates = list(diann_dir.glob("*.protein_description.tsv"))

    if len(candidates) == 0:
        raise FileNotFoundError(
            f"No protein_description.tsv found in: {diann_dir}"
        )

    if len(candidates) == 1:
        return candidates[0]

    raise FileNotFoundError(
        f"Could not uniquely determine protein_description.tsv in: {diann_dir}\n"
        f"Candidates found: {[p.name for p in candidates]}"
    )

# Protein description loading

def load_protein_sequences(protein_desc_file: Path, mutant_prefix: str):
    prot = pd.read_csv(protein_desc_file, sep="\t")

    id_col = find_first_existing(prot, PROTSEQ_ID_COLUMNS)
    name_col = find_first_existing(prot, PROTSEQ_NAME_COLUMNS)
    gene_col = find_first_existing(prot, PROTSEQ_GENE_COLUMNS)
    seq_col = find_first_existing(prot, PROTSEQ_SEQUENCE_COLUMNS)

    if id_col is None:
        raise ValueError(
            f"Could not find a protein ID column in {protein_desc_file}. Tried: {PROTSEQ_ID_COLUMNS}"
        )
    if seq_col is None:
        raise ValueError(
            f"Could not find a sequence column in {protein_desc_file}. Tried: {PROTSEQ_SEQUENCE_COLUMNS}"
        )

    prot = prot.copy()
    prot[seq_col] = prot[seq_col].map(clean_sequence)
    prot = prot[prot[seq_col] != ""].copy()

    prefix_check_cols = [id_col]
    if name_col is not None:
        prefix_check_cols.append(name_col)
    if gene_col is not None:
        prefix_check_cols.append(gene_col)

    prot["is_mutant"] = prot.apply(
        lambda row: row_has_mutant_prefix(row, prefix_check_cols, mutant_prefix),
        axis=1
    )

    mutant_prot = prot[prot["is_mutant"]].copy()
    wt_prot = prot[~prot["is_mutant"]].copy()

    return {
        "all": prot,
        "mutant": mutant_prot,
        "wt": wt_prot,
        "id_col": id_col,
        "name_col": name_col,
        "gene_col": gene_col,
        "seq_col": seq_col,
    }


# Classification

def classify_peptides_against_proteins(peptides, mutant_prot, wt_prot, seq_col, id_col):
    rows = []

    for pep in sorted(set(peptides)):
        pep = clean_sequence(pep)
        if not pep:
            continue

        mutant_hits = find_matching_proteins(pep, mutant_prot, seq_col, id_col)
        wt_hits = find_matching_proteins(pep, wt_prot, seq_col, id_col)

        rows.append({
            "Stripped.Sequence": pep,
            "in_mutant_proteins": len(mutant_hits) > 0,
            "in_wt_proteins": len(wt_hits) > 0,
            "n_mutant_proteins_matched": len(mutant_hits),
            "n_wt_proteins_matched": len(wt_hits),
            "mutant_protein_ids": ";".join(mutant_hits),
            "wt_protein_ids": ";".join(wt_hits),
            "mutation_specific": (len(mutant_hits) > 0 and len(wt_hits) == 0),
        })

    return pd.DataFrame(
        rows,
        columns=[
            "Stripped.Sequence",
            "in_mutant_proteins",
            "in_wt_proteins",
            "n_mutant_proteins_matched",
            "n_wt_proteins_matched",
            "mutant_protein_ids",
            "wt_protein_ids",
            "mutation_specific",
        ]
    )

def summarize_peptides(df: pd.DataFrame, quantity_col: str | None) -> pd.DataFrame:
    if "Stripped.Sequence" not in df.columns:
        raise ValueError("Missing 'Stripped.Sequence' column in report data.")

    agg = {}
    if "Run" in df.columns:
        agg["Run"] = pd.Series.nunique
    if "Precursor.Id" in df.columns:
        agg["Precursor.Id"] = pd.Series.nunique
    if "Modified.Sequence" in df.columns:
        agg["Modified.Sequence"] = lambda x: ";".join(sorted(set(map(str, x.dropna()))))
    if "Protein.Ids" in df.columns:
        agg["Protein.Ids"] = lambda x: ";".join(sorted(set(map(str, x.dropna()))))
    if "Protein.Names" in df.columns:
        agg["Protein.Names"] = lambda x: ";".join(sorted(set(map(str, x.dropna()))))
    if "Genes" in df.columns:
        agg["Genes"] = lambda x: ";".join(sorted(set(map(str, x.dropna()))))

    summary = df.groupby("Stripped.Sequence", dropna=False).agg(agg).reset_index()

    rename_map = {}
    if "Run" in summary.columns:
        rename_map["Run"] = "n_runs"
    if "Precursor.Id" in summary.columns:
        rename_map["Precursor.Id"] = "n_precursors"
    summary = summary.rename(columns=rename_map)

    if quantity_col and quantity_col in df.columns:
        qdf = (
            df.groupby("Stripped.Sequence")[quantity_col]
            .agg(["max", "mean", "median"])
            .reset_index()
            .rename(columns={
                "max": f"{quantity_col}.max",
                "mean": f"{quantity_col}.mean",
                "median": f"{quantity_col}.median",
            })
        )
        summary = summary.merge(qdf, on="Stripped.Sequence", how="left")

    return summary


def main():
    parser = argparse.ArgumentParser(
        description="Find mutation-specific peptides from a DIA-NN output folder."
    )
    parser.add_argument(
        "diann_dir",
        help="Path to DIA-NN output folder"
    )
    parser.add_argument(
        "--mutant-prefix",
        default=DEFAULT_MUTANT_PREFIX,
        help=f"Prefix marking mutant proteins/genes (default: {DEFAULT_MUTANT_PREFIX})"
    )
    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory. Default: <diann_dir>/mutation_specific_results"
    )
    parser.add_argument(
        "--max-qvalue",
        type=float,
        default=DEFAULT_MAX_QVALUE,
        help=f"Max Q.Value threshold (default: {DEFAULT_MAX_QVALUE})"
    )
    parser.add_argument(
        "--max-global-qvalue",
        type=float,
        default=DEFAULT_MAX_GLOBAL_QVALUE,
        help=f"Max Global.Q.Value threshold (default: {DEFAULT_MAX_GLOBAL_QVALUE})"
    )
    args = parser.parse_args()

    diann_dir = Path(args.diann_dir).resolve()
    if not diann_dir.exists() or not diann_dir.is_dir():
        print(f"Error: DIA-NN folder does not exist or is not a directory: {diann_dir}", file=sys.stderr)
        sys.exit(1)

    outdir = Path(args.outdir).resolve() if args.outdir else diann_dir / "mutation_specific_results"
    outdir.mkdir(parents=True, exist_ok=True)

    report_file = discover_report_file(diann_dir)
    protein_desc_file = discover_protein_description_file(diann_dir)

    print(f"DIA-NN folder: {diann_dir}")
    print(f"Main report: {report_file.name}")
    print(f"Protein description: {protein_desc_file.name}")
    print(f"Mutant prefix: {args.mutant_prefix}")
    print(f"Output dir: {outdir}")

    print("\nReading DIA-NN report...")
    report = pd.read_parquet(report_file)
    print(f"report rows: {len(report):,}")
    print(f"report columns: {len(report.columns):,}")

    quantity_col = choose_quantity_column(report)
    print(f"quantity column used: {quantity_col}")

    print("\nFiltering mutant-tagged rows in report.parquet...")
    mutant_mask = build_mutant_tag_mask(report, REPORT_SEARCH_COLUMNS, args.mutant_prefix)
    tagged = report.loc[mutant_mask].copy()
    print(f"rows matching '{args.mutant_prefix}': {len(tagged):,}")

    tagged = filter_qvalues(tagged, args.max_qvalue, args.max_global_qvalue)
    print(f"rows after q-value filtering: {len(tagged):,}")

    keep_cols = existing_columns(tagged, REPORT_KEEP_COLUMNS)
    tagged = tagged[keep_cols].copy()

    tagged_out = outdir / "candidate_mutant_tagged_precursors.tsv"
    tagged.to_csv(tagged_out, sep="\t", index=False)
    print(f"saved: {tagged_out}")

    if "Stripped.Sequence" not in tagged.columns:
        raise ValueError("Column 'Stripped.Sequence' not found in main DIA-NN report")

    print("\nReading protein sequences from protein_description.tsv ...")
    seqinfo = load_protein_sequences(protein_desc_file, args.mutant_prefix)

    mutant_prot = seqinfo["mutant"]
    wt_prot = seqinfo["wt"]
    id_col = seqinfo["id_col"]
    seq_col = seqinfo["seq_col"]

    print(f"protein sequence ID column: {id_col}")
    print(f"protein sequence column: {seq_col}")
    print(f"mutant proteins found: {len(mutant_prot):,}")
    print(f"WT/canonical proteins found: {len(wt_prot):,}")

    mutant_prot.to_csv(outdir / "mutant_proteins_from_description.tsv", sep="\t", index=False)

    print("\nClassifying stripped peptides against mutant/WT protein sequences...")
    peptides = tagged["Stripped.Sequence"].dropna().astype(str).map(clean_sequence)
    peptides = [x for x in peptides if x]

    class_df = classify_peptides_against_proteins(
        peptides=peptides,
        mutant_prot=mutant_prot,
        wt_prot=wt_prot,
        seq_col=seq_col,
        id_col=id_col,
    )

    class_out = outdir / "all_tagged_peptides_classification.tsv"
    class_df.to_csv(class_out, sep="\t", index=False)
    print(f"saved: {class_out}")

    if class_df.empty:
        print("No candidate peptides remained after mutant-tag filtering and sequence cleaning.")
        print("Please check:")
        print("1. whether report.parquet actually contains mutant-tagged rows")
        print("2. whether Stripped.Sequence exists and is non-empty")
        print("3. whether the mutant prefix matches your DIA-NN annotation")
        class_df.to_csv(class_out, sep="\t", index=False)
        return
    ms_peptides = class_df[class_df["mutation_specific"]].copy()
    print(f"mutation-specific peptides found: {len(ms_peptides):,}")

    ms_set = set(ms_peptides["Stripped.Sequence"].tolist())
    ms_precursors = tagged[tagged["Stripped.Sequence"].map(clean_sequence).isin(ms_set)].copy()

    ms_summary = summarize_peptides(ms_precursors, quantity_col=quantity_col)
    ms_summary = ms_summary.merge(ms_peptides, on="Stripped.Sequence", how="left")

    sort_col = None
    if quantity_col and f"{quantity_col}.max" in ms_summary.columns:
        sort_col = f"{quantity_col}.max"
    elif "n_runs" in ms_summary.columns:
        sort_col = "n_runs"

    if sort_col:
        ms_summary = ms_summary.sort_values(sort_col, ascending=False)

    peptide_out = outdir / "mutation_specific_peptides.tsv"
    precursor_out = outdir / "mutation_specific_precursors.tsv"

    ms_summary.to_csv(peptide_out, sep="\t", index=False)
    ms_precursors.to_csv(precursor_out, sep="\t", index=False)

    print(f"saved: {peptide_out}")
    print(f"saved: {precursor_out}")

    print("\nTop mutation-specific peptides:")
    with pd.option_context("display.max_columns", None, "display.width", 220):
        print(ms_summary.head(20))

    print("\nDone.")


if __name__ == "__main__":
    main()
