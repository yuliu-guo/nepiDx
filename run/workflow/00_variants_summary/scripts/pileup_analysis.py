#!/usr/bin/env python3
"""
pileup_analysis.py

Fast allele-support counting for SNVs + simple anchored indels (VCF-style),
with support for opening BAM/CRAM once and reusing a handle.

- open_alignment(): open BAM/CRAM once - pysam.AlignmentFile
- count_allele_support(): SNVs use count_coverage; indels use pileup indel events
- Backward-compatible helpers: get_pileup_using_pysam(), filter_by_quality_score(), etc.
"""

import os
import warnings
from collections import Counter
import pysam

# Suppress a common CRAM-related warning spam (optional)
warnings.filterwarnings("ignore", message="multiple_iterators not implemented for CRAM")

# Reference FASTA needed for CRAM (and ok to provide for BAM too)
REF_FASTA = (
    "/ref/h.sapiens/hg38/genomes/GCA_000001405.15_GRCh38_no_alt_analysis_set_maskedGRC_exclusions.fasta"
)



# Open file once
def open_alignment(bam_path: str, reference_fasta: str = REF_FASTA) -> pysam.AlignmentFile:
    """
    Open a BAM/CRAM and return a pysam.AlignmentFile handle.
    Caller is responsible for closing (bam_handle.close()).
    """
    if not bam_path or not os.path.exists(bam_path):
        raise FileNotFoundError(f"Alignment file not found: {bam_path!r}")
    mode = "rc" if str(bam_path).endswith(".cram") else "rb"
    return pysam.AlignmentFile(bam_path, mode, reference_filename=reference_fasta)



# Classfy variants
def classify_variant(ref: str, alt: str) -> str:
    """
    Simple classifier for allele type.
    - snv: len(ref)==1 and len(alt)==1
    - insertion: len(alt)>len(ref)
    - deletion: len(ref)>len(alt)
    - complex: other equal-length or anything unexpected
    """
    if ref is None or alt is None:
        return "unknown"
    if len(ref) == 1 and len(alt) == 1:
        return "snv"
    if len(alt) > len(ref):
        return "insertion"
    if len(ref) > len(alt):
        return "deletion"
    return "complex"


def count_allele_support(bam_handle: pysam.AlignmentFile,
                         chrom: str,
                         pos: int,
                         ref: str,
                         alt: str,
                         min_bq: int = 30):
    """
    Returns (vaf, alt_n, total_n) for:
      - SNVs (fast path): uses bam_handle.count_coverage
      - simple VCF-anchored insertions/deletions: uses pileup indel events

    VCF anchoring reminder:
      insertion example:  REF=A, ALT=AT, pos=3004  => insertion length 1 after anchor base at 3004
      deletion example:   REF=AGT, ALT=A, pos=3004 => deletion length 2 after anchor base at 3004

    Notes:
    - For indels, alt_n counts reads with the matching indel length at the anchor position.
    - "total_n" is number of reads counted passing base quality at the anchor base (for indels),
      or sum(A,C,G,T) at the SNV position for SNVs.
    - VAF is computed as alt / (alt + ref_like) where ref_like is reads with no indel at that site (indels),
      or alt/(alt+ref) for SNVs.
    """
    vtype = classify_variant(ref, alt)

    # SNV fast path
    if vtype == "snv":
        # count_coverage uses 0-based half-open intervals
        A, C, G, T = bam_handle.count_coverage(chrom, pos - 1, pos, quality_threshold=min_bq)
        counts = {"A": A[0], "C": C[0], "G": G[0], "T": T[0]}
        total = sum(counts.values())
        ref_n = counts.get(ref, 0)
        alt_n = counts.get(alt, 0)
        informative = ref_n + alt_n
        vaf = (alt_n / informative) if informative else 0.0
        return vaf, alt_n, total

    # Indels: use pileup indel events at anchor base
    ins_len = 0
    del_len = 0
    if vtype == "insertion":
        ins_len = len(alt) - len(ref)
    elif vtype == "deletion":
        del_len = len(ref) - len(alt)
    else:
        
        # return 0 alt support but still compute total base coverage at anchor.
        A, C, G, T = bam_handle.count_coverage(chrom, pos - 1, pos, quality_threshold=min_bq)
        total = int(A[0] + C[0] + G[0] + T[0])
        return 0.0, 0, total

    total = 0
    alt_n = 0
    ref_like_n = 0  # reads with no indel at the anchor

    for col in bam_handle.pileup(chrom, pos - 1, pos, truncate=True):
        if col.pos != pos - 1:
            continue
        for pr in col.pileups:
            if pr.is_refskip:
                continue

            # use anchor base quality filtering; if anchor is deleted in the read at this column,
            # query_position is None. We skip those reads to avoid inflating totals.
            if pr.is_del:
                continue

            qpos = pr.query_position
            aln = pr.alignment
            bq = aln.query_qualities[qpos]
            if bq < min_bq:
                continue

            total += 1

            # pr.indel: +k insertion after this base, -k deletion after this base, 0 no indel
            if ins_len and pr.indel == ins_len:
                alt_n += 1
            elif del_len and pr.indel == -del_len:
                alt_n += 1
            elif pr.indel == 0:
                ref_like_n += 1

    informative = alt_n + ref_like_n
    vaf = (alt_n / informative) if informative else 0.0
    return vaf, alt_n, total


def pileup_counts_and_vaf(bam_or_handle,
                          chrom: str,
                          location: int,
                          ref: str,
                          alt: str,
                          quality_threshold: int = 30,
                          reference_fasta: str = REF_FASTA):
    """
    Wrapper that accepts either a BAM path or an open handle and returns:
      (vaf, alt_n, total_n)

    """
    close_after = False 
    if isinstance(bam_or_handle, str):
        samfile = open_alignment(bam_or_handle, reference_fasta=reference_fasta)
        close_after = True
    else:
        samfile = bam_or_handle

    vaf, alt_n, total_n = count_allele_support(
        samfile, chrom, location, ref, alt, min_bq=quality_threshold
    )

    if close_after:
        samfile.close()

    return vaf, alt_n, total_n
    
    

