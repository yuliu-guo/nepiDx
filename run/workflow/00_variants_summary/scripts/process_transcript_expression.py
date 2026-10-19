import pandas as pd

# file_path -> dict(Name -> (TPM, NumReads))
_QUANT_LOOKUP_CACHE = {}

def load_quant_lookup(file_path):
    """Load quant.sf once into a lookup dict."""
    if file_path in _QUANT_LOOKUP_CACHE:
        return _QUANT_LOOKUP_CACHE[file_path]

    df = pd.read_csv(file_path, sep="\t", usecols=["Name", "TPM", "NumReads"])
    lookup = {row["Name"]: (float(row["TPM"]), float(row["NumReads"])) for _, row in df.iterrows()}

    _QUANT_LOOKUP_CACHE[file_path] = lookup
    return lookup


def get_tpm(file_path, refseq_id):
    lookup = load_quant_lookup(file_path)
    v = lookup.get(refseq_id)
    return v[0] if v else None


def get_numreads(file_path, refseq_id):
    lookup = load_quant_lookup(file_path)
    v = lookup.get(refseq_id)
    return v[1] if v else None