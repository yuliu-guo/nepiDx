#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
from query_best_match import BestMatcher


# count runtime
import time
from collections import defaultdict

# This script only works for SV neoepitope summary 

timings = defaultdict(float)

def _tstart():
    return time.perf_counter()

def _tadd(key, t0):
    timings[key] += (time.perf_counter() - t0)

def split_kmers(cell_value: str, sep: str):
    """
    seq-col already contains kmers separated by sep (default ';').
    Returns a list of kmers (uppercased, stripped), skipping empties.
    """
    if cell_value is None:
        return []
    s = str(cell_value).strip()
    if not s or s.lower() == "nan":
        return []
    parts = [p.strip().upper() for p in s.split(sep)]
    return [p for p in parts if p]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-tsv", required=True, help="Input TSV path")
    ap.add_argument("--output-tsv", required=True, help="Output TSV path")

    # flexible input columns
    ap.add_argument("--seq-col", required=True, help="Column containing pre-sliced kmers separated by a delimiter")
    ap.add_argument(
        "--keep-cols",
        nargs="*",
        default=[],
        help="Column names to copy into each expanded output row"
    )

    ap.add_argument("--kmer", type=int, default=9, help="k-mer size (default: 9) (used only for filtering)")
    ap.add_argument("--sep", default=";", help="Delimiter for kmers in --seq-col (default: ';')")

    ap.add_argument("--db", required=True, help="Path to human_k*.sqlite")
    ap.add_argument("--blos", required=True, help="Path to BLOSUM62-2 file")
    ap.add_argument("--beta", type=float, default=0.125, help="beta for kernel")
    ap.add_argument("--top-n-seed-hits", type=int, default=200, help="Candidate cap for WT search")
    ap.add_argument("--n-finalists", type=int, default=20, help="Kernel is only computed on these top Hamming hits")
    ap.add_argument("--cache-path", default=None, help="Optional K1 cache path")
    ap.add_argument("--p-min", type=int, default=None, help="Minimum substring length for K3 (default: kmer)")
    ap.add_argument("--p-max", type=int, default=None, help="Maximum substring length for K3 (default: kmer)")
    args = ap.parse_args()

    k = args.kmer
    p_min = args.p_min if args.p_min is not None else k
    p_max = args.p_max if args.p_max is not None else k

    t0 = _tstart()
    df = pd.read_csv(args.input_tsv, sep="\t", dtype=str)
    _tadd("read_csv", t0)

    # validate columns
    if args.seq_col not in df.columns:
        raise ValueError(f"--seq-col '{args.seq_col}' not found in input TSV")

    missing_keep = [c for c in args.keep_cols if c not in df.columns]
    if missing_keep:
        raise ValueError(f"--keep-cols not found in input TSV: {missing_keep}")

    out_rows = []

    t0 = _tstart()
    matcher = BestMatcher(
        db_path=args.db,
        blos_path=args.blos,
        beta=args.beta,
        cache_path=args.cache_path,
        p_min=p_min,
        p_max=p_max,
        n_finalists=args.n_finalists,
    )
    _tadd("init_matcher_once", t0)

    try:
        for _, row in df.iterrows():
            cell = row[args.seq_col]

            t_split = _tstart()
            kmers = split_kmers(cell, args.sep)
            _tadd("split_kmers", t_split)

            # Skip rows with no kmers
            if not kmers:
                _tadd("skipped_empty_rows", _tstart())  # tiny marker
                continue

            # copy requested columns
            base_info = {col: row[col] for col in args.keep_cols}

            t_all = _tstart()
            for mut_peptide in kmers:
                # optional length filter
                if k and len(mut_peptide) != k:
                    continue

                t_one = _tstart()
                result = matcher.find_best_match(
                    peptide=mut_peptide,
                    top_n_seed_hits=args.top_n_seed_hits,
                )

                out_row = dict(base_info)
                out_row[args.seq_col] = cell  # keep original cell content as-is
                out_row["Mut_peptide"] = mut_peptide
                out_row["WT_peptide"] = result["wt_peptide"]
                out_row["hamming_similarity"] = result["hamming_similarity"]
                out_row["blosum62_2_similarity"] = result["blosum62_2_similarity"]

                out_rows.append(out_row)
                _tadd("find_best_match_one", t_one)

            _tadd("find_best_match_all", t_all)

    finally:
        t0 = _tstart()
        matcher.close()
        _tadd("close_matcher_once", t0)


    ordered_cols = (
        list(args.keep_cols)
        + [args.seq_col, "Mut_peptide", "WT_peptide", "hamming_similarity", "blosum62_2_similarity"]
    )
    
    if not out_rows:
        print(
            "WARNING: There are no available kmers; writing empty output TSV with headers only.",
            file=sys.stderr
        )
        
    out_df = pd.DataFrame(out_rows, columns=ordered_cols)

    t0 = _tstart()
    out_df.to_csv(args.output_tsv, sep="\t", index=False)
    _tadd("write_tsv", t0)

    print("=== Timing breakdown (seconds) ===")
    for key, val in sorted(timings.items(), key=lambda kv: -kv[1]):
        print(f"{key:28s} {val:10.3f}")


if __name__ == "__main__":
    main()
