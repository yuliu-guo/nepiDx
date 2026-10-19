#!/usr/bin/env python3
import argparse
import json
import sqlite3

import hamming_similarity
from blosum62_2_kernel_similarity import KernelScorer



# SQLite


def get_meta(conn, key):
    cur = conn.cursor()
    row = cur.execute("SELECT value FROM meta WHERE key = ?", (key,)).fetchone()
    if row is None:
        raise ValueError(f"Missing meta key: {key}")
    return row[0]


def seeds_of(seq, seed_len):
    seen = set()
    for i in range(len(seq) - seed_len + 1):
        s = seq[i:i + seed_len]
        if s not in seen:
            seen.add(s)
            yield s


def fetch_candidates(conn, peptide, seed_len, top_n_seed_hits=200):
    cur = conn.cursor()
    seed_list = list(seeds_of(peptide, seed_len))
    if not seed_list:
        return []

    placeholders = ",".join("?" for _ in seed_list)
    sql = f"""
        SELECT k.seq, COUNT(*) AS shared_seed_count
        FROM seeds s
        JOIN kmers k ON k.id = s.kmer_id
        WHERE s.seed IN ({placeholders})
        GROUP BY k.id
        ORDER BY shared_seed_count DESC
        LIMIT ?
    """
    rows = cur.execute(sql, (*seed_list, top_n_seed_hits)).fetchall()
    return [r[0] for r in rows]


def best_match_two_stage(
    conn,
    peptide,
    seed_len,
    kernel_scorer,
    top_n_seed_hits=200,
    n_finalists=5,
    p_min=9,
    p_max=9,
):
    peptide = peptide.strip().upper()

    candidates = fetch_candidates(conn, peptide, seed_len, top_n_seed_hits=top_n_seed_hits)

    if not candidates:
        cur = conn.cursor()
        candidates = [r[0] for r in cur.execute("SELECT seq FROM kmers")]

    if not candidates:
        return {
            "wt_peptide": None,
            "hamming_similarity": None,   # will store Hamming similarity
            "blosum62_2_similarity": None,
        }

    #  compute Hamming for all candidates
    ham_scored = []
    for wt in candidates:
        ham = hamming_similarity.hamming_similarity(peptide, wt)
        ham_scored.append((wt, ham))

    # Find the 5th-highest Hamming score 
    ham_values = sorted((h for _, h in ham_scored), reverse=True)
    cutoff = ham_values[min(n_finalists - 1, len(ham_values) - 1)]

    # Keep ALL peptides with ham >= cutoff 
    finalists = [(wt, ham) for wt, ham in ham_scored if ham >= cutoff]

    # only compute kernel for the highest top5 hamming score
    best = None
    best_key = None

    for wt, ham in finalists:
        ker = kernel_scorer.k3_hat(peptide, wt, p_min=p_min, p_max=p_max)
        key = (-ham, -ker, wt)
        if best is None or key < best_key:
            best = {
                "wt_peptide": wt,
                "hamming_similarity": float(ham),
                "blosum62_2_similarity": float(ker),
            }
            best_key = key

    return best

class BestMatcher:
    """
    Keep the DB connection and KernelScorer alive across many queries.
      1) Hamming first
      2) kernel only on finalists
    """

    def __init__(
        self,
        db_path,
        blos_path,
        beta=0.125,
        cache_path=None,
        p_min=9,
        p_max=9,
        n_finalists=5,
    ):
        self.conn = sqlite3.connect(db_path)
        self.seed_len = int(get_meta(self.conn, "seed_len"))
        self.k = int(get_meta(self.conn, "k"))

        self.kernel_scorer = KernelScorer(
            blos_path=blos_path,
            beta=beta,
            cache_path=cache_path,
        )

        self.p_min = p_min
        self.p_max = p_max
        self.n_finalists = n_finalists

    def close(self):
        self.conn.close()

    def find_best_match(self, peptide, top_n_seed_hits=200):
        peptide = peptide.strip().upper()

        if len(peptide) != self.k:
            raise ValueError(
                f"Query peptide length is {len(peptide)}, but DB was built with k={self.k}"
            )

        return best_match_two_stage(
            conn=self.conn,
            peptide=peptide,
            seed_len=self.seed_len,
            kernel_scorer=self.kernel_scorer,
            top_n_seed_hits=top_n_seed_hits,
            n_finalists=self.n_finalists,
            p_min=self.p_min,
            p_max=self.p_max,
        )



def find_best_match(
    db_path,
    blos_path,
    peptide,
    beta=0.125,
    top_n_seed_hits=200,
    cache_path=None,
    p_min=9,
    p_max=9,
    n_finalists=5,
):
    matcher = BestMatcher(
        db_path=db_path,
        blos_path=blos_path,
        beta=beta,
        cache_path=cache_path,
        p_min=p_min,
        p_max=p_max,
        n_finalists=n_finalists,
    )
    try:
        return matcher.find_best_match(
            peptide=peptide,
            top_n_seed_hits=top_n_seed_hits,
        )
    finally:
        matcher.close()



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", required=True, help="SQLite DB built from WT k-mers")
    ap.add_argument("--blos", required=True, help="Path to BLOSUM62-2 file")
    ap.add_argument("--peptide", required=True, help="Query peptide")
    ap.add_argument("--beta", type=float, default=0.125, help="beta for K1=(BLOSUM62-2)^beta")
    ap.add_argument("--top-n-seed-hits", type=int, default=200, help="Candidate cap before Hamming stage")
    ap.add_argument("--n-finalists", type=int, default=20, help="How many top Hamming hits get kernel scoring")
    ap.add_argument("--cache-path", default=None, help="Optional explicit path for cached K1")
    ap.add_argument("--p-min", type=int, default=9, help="Minimum substring length for K3")
    ap.add_argument("--p-max", type=int, default=9, help="Maximum substring length for K3")
    args = ap.parse_args()

    matcher = BestMatcher(
        db_path=args.db,
        blos_path=args.blos,
        beta=args.beta,
        cache_path=args.cache_path,
        p_min=args.p_min,
        p_max=args.p_max,
        n_finalists=args.n_finalists,
    )
    try:
        result = matcher.find_best_match(
            peptide=args.peptide,
            top_n_seed_hits=args.top_n_seed_hits,
        )
    finally:
        matcher.close()

    print(json.dumps(result, ensure_ascii=False))


if __name__ == "__main__":
    main()