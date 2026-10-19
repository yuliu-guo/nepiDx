import os


# Reference file paths
PROTEIN_FASTA = "/ref/GRCh38_latest_protein.faa"
CCDS_REFSEQ_TABLE = "/ref/CCDS2Sequence.current.txt"
CCDS_NUCLEOTIDE_FASTA = "/ref/CCDS_nucleotide.current.fna"


# Module-level lookups 
_PROTEIN_LOOKUP = None          # dict: base_ENSP -> protein sequence
_CCDS_ID_LOOKUP = None          # dict: base_RefSeq -> CCDS_id
_CCDS_SEQ_LOOKUP = None         # dict: base_CCDS_id -> nucleotide sequence


def _base_id(x: str) -> str:
    """Return ID without version suffix. E.g. NP_XXXXX.2 -> NP_XXXXX"""
    return str(x).split(".", 1)[0]



# ENSP to sequence
def _build_protein_lookup(fasta_path: str) -> dict:
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"Protein FASTA not found: {fasta_path}")

    lookup = {}
    seq_id = None
    seq_chunks = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # flush previous record
                if seq_id is not None:
                    base = _base_id(seq_id)
                    # keep first occurrence for that base id
                    lookup.setdefault(base, "".join(seq_chunks))

                # start new record
                seq_id = line[1:].split()[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        # flush last record
        if seq_id is not None:
            base = _base_id(seq_id)
            lookup.setdefault(base, "".join(seq_chunks))

    return lookup


def get_protein_seq(ensp_id):
    """Return protein sequence for ENSP (version-agnostic)."""
    if not ensp_id:
        return None

    global _PROTEIN_LOOKUP
    if _PROTEIN_LOOKUP is None:
        _PROTEIN_LOOKUP = _build_protein_lookup(PROTEIN_FASTA)

    return _PROTEIN_LOOKUP.get(_base_id(ensp_id))



#  RefSeq transcript to CCDS ID
def _build_ccds_id_lookup(table_path: str) -> dict:
    if not os.path.exists(table_path):
        raise FileNotFoundError(f"CCDS2Sequence table not found: {table_path}")

    lookup = {}
    with open(table_path, "r") as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            # RefSeq transcript is column 5 (index 4), CCDS is column 1 (index 0)
            if len(parts) <= 4:
                continue
            ccds_id = parts[0]
            refseq = parts[4]
            if not refseq:
                continue
            lookup.setdefault(_base_id(refseq), ccds_id)
    return lookup


def get_CCDS_ID(RefSeq_transcript_id):
    """Return CCDS ID for a RefSeq transcript ID (version-agnostic)."""
    if not RefSeq_transcript_id:
        return None

    global _CCDS_ID_LOOKUP
    if _CCDS_ID_LOOKUP is None:
        _CCDS_ID_LOOKUP = _build_ccds_id_lookup(CCDS_REFSEQ_TABLE)

    return _CCDS_ID_LOOKUP.get(_base_id(RefSeq_transcript_id))


# CCDS ID to nucleotide sequence

def _build_ccds_seq_lookup(fasta_path: str) -> dict:
    if not os.path.exists(fasta_path):
        raise FileNotFoundError(f"CCDS nucleotide FASTA not found: {fasta_path}")

    lookup = {}
    seq_id = None
    seq_chunks = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                # flush previous record
                if seq_id is not None:
                    base = _base_id(seq_id)
                    lookup.setdefault(base, "".join(seq_chunks))

                # header format: >CCDSxxxx|...
                seq_id = line[1:].split("|")[0]
                seq_chunks = []
            else:
                seq_chunks.append(line)

        # flush last record
        if seq_id is not None:
            base = _base_id(seq_id)
            lookup.setdefault(base, "".join(seq_chunks))

    return lookup


def get_CCDS_seq(CCDS_id):
    """Return CCDS nucleotide sequence for CCDS ID (version-agnostic)."""
    if not CCDS_id:
        return None

    global _CCDS_SEQ_LOOKUP
    if _CCDS_SEQ_LOOKUP is None:
        _CCDS_SEQ_LOOKUP = _build_ccds_seq_lookup(CCDS_NUCLEOTIDE_FASTA)

    return _CCDS_SEQ_LOOKUP.get(_base_id(CCDS_id))
    
    
