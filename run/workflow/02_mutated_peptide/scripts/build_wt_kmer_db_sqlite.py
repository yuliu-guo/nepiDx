import argparse
import sqlite3
from pathlib import Path


def read_fasta(path):
    name = None
    seq_chunks = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq_chunks)
                name = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line.upper())
        if name is not None:
            yield name, "".join(seq_chunks)


def kmers(seq, k):
    n = len(seq)
    for i in range(n - k + 1):
        mer = seq[i:i+k]
        if mer.isalpha():
            yield mer


def seeds_of(seq, seed_len):
    # unique seeds from one peptide
    seen = set()
    for i in range(len(seq) - seed_len + 1):
        
        s = seq[i:i+seed_len]
        if s not in seen:
            seen.add(s)
            yield s


def init_db(conn):
    cur = conn.cursor()
    cur.execute("""
        CREATE TABLE IF NOT EXISTS meta (
            key TEXT PRIMARY KEY,
            value TEXT
        )
    """)
    cur.execute("""
        CREATE TABLE IF NOT EXISTS kmers (
            id INTEGER PRIMARY KEY,
            seq TEXT UNIQUE
        )
    """)
    cur.execute("""
        CREATE TABLE IF NOT EXISTS seeds (
            seed TEXT,
            kmer_id INTEGER,
            UNIQUE(seed, kmer_id)
        )
    """)
    cur.execute("CREATE INDEX IF NOT EXISTS idx_kmers_seq ON kmers(seq)")
    cur.execute("CREATE INDEX IF NOT EXISTS idx_seeds_seed ON seeds(seed)")
    conn.commit()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fasta", required=True, help="WT proteome FASTA")
    ap.add_argument("--db", required=True, help="Output SQLite DB")
    ap.add_argument("--k", type=int, default=9, help="WT k-mer size")
    ap.add_argument("--seed-len", type=int, default=4, help="Seed size for candidate lookup")
    args = ap.parse_args()

    if args.seed_len > args.k:
        raise ValueError("--seed-len cannot be bigger than --k")

    db_path = Path(args.db)
    conn = sqlite3.connect(db_path)
    init_db(conn)
    cur = conn.cursor()

    cur.execute("DELETE FROM meta")
    cur.execute("INSERT INTO meta(key, value) VALUES (?, ?)", ("k", str(args.k)))
    cur.execute("INSERT INTO meta(key, value) VALUES (?, ?)", ("seed_len", str(args.seed_len)))
    conn.commit()

    inserted = 0
    for _, seq in read_fasta(args.fasta):
        for mer in kmers(seq, args.k):
            cur.execute("INSERT OR IGNORE INTO kmers(seq) VALUES (?)", (mer,))
            # sqlite rowcount is 1 only if newly inserted
            if cur.rowcount == 1:
                kmer_id = cur.lastrowid
                for sd in seeds_of(mer, args.seed_len):
                    cur.execute(
                        "INSERT OR IGNORE INTO seeds(seed, kmer_id) VALUES (?, ?)",
                        (sd, kmer_id)
                    )
                inserted += 1

    conn.commit()
    print(f"Unique WT {args.k}-mers stored: {inserted}")
    conn.close()


if __name__ == "__main__":
    main()


# example command: python build_wt_kmer_db_sqlite.py --fasta /path/to/protein_fasta/UP000005640_9606.fasta --db human_k9.sqlite --k 9 --seed-len 4