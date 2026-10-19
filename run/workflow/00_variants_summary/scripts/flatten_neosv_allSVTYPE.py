#!/usr/bin/env python3

# Requirement: module load Conda-BioModules/1.1.0 (The core is to installpyensembl)

import os
import re
import argparse
import csv
from warnings import warn
from collections import Counter

import vcfpy
from pyensembl.genome import Genome
#from pyensembl import EnsemblRelease
from Bio.Seq import Seq



# VCF parsing helpers 


def vcf_alt_format_check(alt: str) -> bool:
    """
    True if ALT looks like a breakend string:
      right: <seq> ] chr:pos ]  OR <seq> [ chr:pos [
      left : ] chr:pos ] <seq>  OR [ chr:pos [ <seq>
    """
    legal_pattern_right = re.compile(r'([ATCGatcg]+)(\]|\[)(.+?):(\d+)(\]|\[)$')
    legal_pattern_left = re.compile(r'(\]|\[)(.+?):(\d+)(\]|\[)([ATCGatcg]+)$')
    return bool(legal_pattern_right.match(alt) or legal_pattern_left.match(alt))


def _alt_to_string_vcfpy(alt_obj) -> str:
    """vcfpy ALT objects (BreakEnd, symbolic ALT, etc.) usually have serialize()."""
    if hasattr(alt_obj, "serialize"):
        return alt_obj.serialize()
    return str(alt_obj)


def infer_pattern_from_alt(alt: str) -> int:
    # t[p[
    pattern1 = re.compile(r"[ATCGatcg]+\[.+?\[")
    # t]p]
    pattern2 = re.compile(r"[ATCGatcg]+\].+?\]")
    # ]p]t
    pattern3 = re.compile(r"\].+?\][ATCGatcg]+")
    if pattern1.search(alt):
        return 1
    if pattern2.search(alt):
        return 2
    if pattern3.search(alt):
        return 3
    return 4  # [p[t


def get_coordinate(alt: str):
    coord_pattern = re.compile(r"(\]|\[)(.+?):(\d+)(\]|\[)")
    coord = coord_pattern.findall(alt)[0]
    chrom, pos = coord[1].replace("chr", ""), int(coord[2])
    return chrom, pos


def get_insertion(alt: str):
    """
    Extract inserted sequence from BND ALT if present.
    Example: GTAGTG]chr1:200333897] -> insertion "TAGTG" (drops anchor base).
    """
    nt_pattern = re.compile(r"[ATCGatcg]+")
    right_pattern = re.compile(r"([ATCGatcg]+)(\]|\[)")
    nt = nt_pattern.findall(alt)[0].upper()
    if len(nt) == 1:
        return ""
    if right_pattern.search(alt):
        return nt[1:]
    return nt[:-1]

def normalize_insseq(x) -> str:
    """Convert INFO SVINSSEQ-like values into a clean DNA string."""
    if x is None:
        return ""
    if isinstance(x, (list, tuple)):
        # common vcfpy case: ["A"] or ["TAGTG"]
        x = "".join([str(i) for i in x])
    else:
        x = str(x)

    # keep only DNA letters, uppercase
    x = re.sub(r"[^ACGTacgt]", "", x).upper()
    return x

def sorted_coord(chrom1, pos1, chrom2, pos2):
    bp1 = f"{chrom1}_{pos1}"
    bp2 = f"{chrom2}_{pos2}"
    return tuple(sorted([bp1, bp2]))


def svtype_from_fields(chrom1, pos1, chrom2, pos2, pattern: int) -> str:
    # same mapping as your StructuralVariant.svtype
    if chrom1 != chrom2:
        return "TRA"
    if pos1 < pos2:
        if pattern == 1:
            return "DEL"
        if pattern == 2:
            return "h2hINV"
        if pattern == 3:
            return "DUP"
        return "t2tINV"
    else:
        if pattern == 1:
            return "DUP"
        if pattern == 2:
            return "h2hINV"
        if pattern == 3:
            return "DEL"
        return "t2tINV"


def standardize_sv_from_vcfpy_record(rec, quiet_skip=True):
    """
    Convert one vcfpy record into:
      - list of SV dicts (for BND with multiple ALT)
      - or one SV dict
      - or None (skip)

    STRICT INS policy:
      - keep INS only if INFO contains SVINSSEQ
      - if only LEFT_SVINSSEQ/RIGHT_SVINSSEQ exist -> skip
    """
    svtype = rec.INFO.get("SVTYPE", None)

    #FILTER = rec.FILTER
    svlen_raw = rec.INFO.get("SVLEN")
    # handle missing SVLEN
    if svlen_raw is None:
        svlen = None
    else:
        # SVLEN might be a single int, or a list/tuple like (-1187,) / [-1187]
        if isinstance(svlen_raw, (list, tuple)):
            svlen = abs(int(svlen_raw[0])) if len(svlen_raw) else None
        else:
            svlen = abs(int(svlen_raw))

    
    
    SOMATICSCORE = rec.INFO.get("SOMATICSCORE", None)

    chrom = str(rec.CHROM).replace("chr", "")
    pos1 = int(rec.POS)

    # BND
    if svtype == "BND":
        out = []
        for alt_obj in rec.ALT:
            alt = _alt_to_string_vcfpy(alt_obj).replace("chr", "")
            if not vcf_alt_format_check(alt):
                continue

            chrom2, pos2 = get_coordinate(alt)
            insertion = get_insertion(alt)  # handles insertion embedded in ALT
            pattern = infer_pattern_from_alt(alt)

            sv = {
                "chrom1": chrom, "pos1": pos1,
                "chrom2": chrom2, "pos2": pos2,
                "insertion": insertion,
                "pattern": pattern,
            }
            #sv['FILTER'] = FILTER
            sv["svtype"] = svtype_from_fields(sv["chrom1"], sv["pos1"], sv["chrom2"], sv["pos2"], sv["pattern"])
            sv['SVTYPE'] = svtype
            sv['SVLEN'] = svlen
            sv['SOMATICSCORE'] = SOMATICSCORE
            sv["sorted_coord"] = sorted_coord(sv["chrom1"], sv["pos1"], sv["chrom2"], sv["pos2"])
            out.append(sv)
        return out if out else None

    # DEL / DUP:TANDEM 
    if svtype in ("DEL", "DUP"):
        end = rec.INFO.get("END", None)
        if end is None:
            return None
        pos2 = int(end)

        # In Manta, duplications for SVTYPE=DUP is DUP:TAMDEM where ALT contains TANDEM
        alt_texts = [_alt_to_string_vcfpy(a) for a in rec.ALT]
        is_tandem = any("TANDEM" in t for t in alt_texts)
        if svtype == "DUP" and (not is_tandem):
            if not quiet_skip:
                warn(f"Skipping non-tandem DUP at {rec.CHROM}:{rec.POS}")
            return None

        # Map to the existing pattern convention:
        # pos1<pos2: pattern 1 behaves like DEL; pattern 3 behaves like DUP
        if svtype == "DEL":
            pattern = 1 if pos1 < pos2 else 3
        else:  # DUP:TANDEM
            pattern = 3 if pos1 < pos2 else 1

        insertion = normalize_insseq(rec.INFO.get("SVINSSEQ", ""))

        sv = {
            "chrom1": chrom, "pos1": pos1,
            "chrom2": chrom, "pos2": pos2,
            "insertion": insertion,
            "pattern": pattern,
        }
        #sv['FILTER'] = FILTER
        sv["svtype"] = svtype_from_fields(sv["chrom1"], sv["pos1"], sv["chrom2"], sv["pos2"], sv["pattern"])
        sv['SVTYPE'] = svtype
        sv['SVLEN'] = svlen
        sv['SOMATICSCORE'] = SOMATICSCORE
        sv["sorted_coord"] = sorted_coord(sv["chrom1"], sv["pos1"], sv["chrom2"], sv["pos2"])
        return sv

    # INS
    if svtype == "INS":
        insseq = rec.INFO.get("SVINSSEQ", None)
        if not insseq:
            # STRICT: skip when only LEFT_SVINSSEQ/RIGHT_SVINSSEQ are present
            # Large insertions are reported in some cases even when the insert 
            # sequence cannot be fully assembled. 
            # In this case, the SV Caller reports the insertion using the <INS>
            # symbolic allele and includes the special INFO fields LEFT_SVINSSEQ and 
            # RIGHT_SVINSSEQ to describe the assembled left and right ends of the insert sequence. 

            return None

        end = rec.INFO.get("END", pos1)
        pos2 = int(end)

        sv = {
            "chrom1": chrom, "pos1": pos1,
            "chrom2": chrom, "pos2": pos2,
            "insertion": normalize_insseq(insseq),
            "pattern": 1,  # treat like 5' + INS + 3'
            "svtype": "INS",
        }
        
        sv['SVTYPE'] = svtype
        sv['SVLEN'] = svlen
        sv['SOMATICSCORE'] = SOMATICSCORE
        sv["sorted_coord"] = sorted_coord(sv["chrom1"], sv["pos1"], sv["chrom2"], sv["pos2"])

        return sv

    return None


def vcf_load(vcf_path: str, quiet_skip: bool = True):
    """
    Load SVs from VCF/VCF.GZ and return standardized SV dicts:
      {"chrom1","pos1","chrom2","pos2","insertion","pattern","svtype","SVTYPE","SVLEN","FILTER","SOMATICSCORE","sorted_coord"}

    Handles BND + DEL + DUP:TANDEM + INS (strict SVINSSEQ-only).
    """
    svs = []
    reader = vcfpy.Reader.from_path(vcf_path)
    print(f'reader is {reader}')

    for rec in reader:
        #print(f'rec is {rec}')
        if 'PASS' not in rec.FILTER:
            #print(f'rec.FILTER is {rec.FILTER}')
            continue
        sv = standardize_sv_from_vcfpy_record(rec, quiet_skip=quiet_skip)
        
        if sv is None:
            continue
        
        if isinstance(sv, list):
            svs.extend(sv)
        else:
            svs.append(sv)

    reader.close()
    print("SV counts:", Counter([s["svtype"] for s in svs]))
    #print("SVS",svs)
    return svs


def remove_duplicate(sv_list):
    seen = set()
    out = []
    for sv in sv_list:
        key = sv["sorted_coord"]
        if key in seen:
            continue
        seen.add(key)
        out.append(sv)
    print(f"{len(sv_list) - len(out)} duplicated SVs were removed")
    return out



# Ensembl + transcript 

def ensembl_load(release: str = "custom"):
    cache_dir = "/ref/pyensembl_cache"
    gtf_file = "/ref/pyensembl_cache/Homo_sapiens.GRCh38.115.gtf.gz"
    cdna_file = "/ref/pyensembl_cache/Homo_sapiens.GRCh38.cdna.all.fa.gz"
    protein_file="/ref/pyensembl_cache/Homo_sapiens.GRCh38.pep.all.fa.gz"
    release = "custom"
    
    ensembl = Genome(
        gtf_path_or_url=gtf_file,
        transcript_fasta_paths_or_urls=[cdna_file],
        protein_fasta_paths_or_urls=[protein_file],
        reference_name="GRCh38.115",
        annotation_name="Ensembl",
    )
    ensembl.index()
    return ensembl


def get_longest_transcript(transcripts):
    transcripts = sorted(transcripts, key=lambda t: t.end - t.start)
    return transcripts[-1]


def get_transcript(chrom, pos, ensembl, complete=True):
    transcripts = ensembl.transcripts_at_locus(contig=str(chrom), position=int(pos))
    #print(f'transcripts is {transcripts}')
    transcripts_comp = [t for t in transcripts if t.complete]
    if transcripts_comp:
        return get_longest_transcript(transcripts_comp)
    if complete:
        return None
    return get_longest_transcript(transcripts) if transcripts else None


def get_cds_range(transcript):
    #print(f'transcript is {transcript}')
    #print(f'coding_sequence_position_ranges is {transcript.coding_sequence_position_ranges}')
    return transcript.coding_sequence_position_ranges


def get_noncds_range(transcript):
    cds_ranges = transcript.coding_sequence_position_ranges
    ncds_ranges = []
    if not cds_ranges:
        return [(0, 0)]

    if transcript.strand == "+":
        utr5_start = transcript.start
        utr5_end = cds_ranges[0][0] - 1
        ncds_ranges.append((0, 0) if utr5_start > utr5_end else (utr5_start, utr5_end))

        for i in range(1, len(cds_ranges)):
            cds_start = cds_ranges[i - 1][1] + 1
            cds_end = cds_ranges[i][0] - 1
            ncds_ranges.append((cds_start, cds_end))

        utr3_start = cds_ranges[-1][1] + 1
        utr3_end = transcript.end
        ncds_ranges.append((0, 0) if utr3_start > utr3_end else (utr3_start, utr3_end))
    else:
        utr5_start = cds_ranges[0][1] + 1
        utr5_end = transcript.end
        ncds_ranges.append((0, 0) if utr5_start > utr5_end else (utr5_start, utr5_end))

        for i in range(1, len(cds_ranges)):
            cds_start = cds_ranges[i][1] + 1
            cds_end = cds_ranges[i - 1][0] - 1
            ncds_ranges.append((cds_start, cds_end))

        utr3_start = transcript.start
        utr3_end = cds_ranges[-1][0] - 1
        ncds_ranges.append((0, 0) if utr3_start > utr3_end else (utr3_start, utr3_end))

    return ncds_ranges


def is_overlap(pos, intervals):
    for i in range(len(intervals)):
        if intervals[i][0] <= pos <= intervals[i][1]:
            return True, i
    return False, -1

def get_gene_name(chrom, pos, ensembl, prefer_biotype="protein_coding") -> str:
    """
    Return gene name(s) at a locus.
    If multiple gene names overlap, returns a ';'-joined deterministic string.
    If none, returns "".
    """
    
    txs = ensembl.transcripts_at_locus(contig=str(chrom), position=int(pos))
    if not txs:
        return ""
    
    # Prefer a gene with the requested biotype if any transcript at locus has it
    for t in txs:
        if (t.biotype or "").lower() == prefer_biotype.lower():
            return t.gene_name

    # Fallback: just take the first gene seen
    return txs[0].gene_name



# mutated-transcript building
def truncate_cds_simple(transcript, direction: str, pos: int):
    cds_ranges = get_cds_range(transcript)
    ncds_ranges = get_noncds_range(transcript)
    if not cds_ranges:
        return None

    cds_overlap, cds_index = is_overlap(pos, cds_ranges)
    ncds_overlap, ncds_index = is_overlap(pos, ncds_ranges)

    # safety: avoid -1 slicing guessing
    if (not cds_overlap) and (not ncds_overlap):
        return None

    selected = []

    if transcript.strand == "+" and direction == "5":
        if cds_overlap:
            selected = list(cds_ranges[: cds_index + 1])
            s, e = selected[-1]
            selected[-1] = (s, pos)
        else:
            selected = list(cds_ranges[: ncds_index])

    elif transcript.strand == "+" and direction == "3":
        if cds_overlap:
            selected = list(cds_ranges[cds_index:])
            s, e = selected[0]
            selected[0] = (pos, e)
        else:
            selected = list(cds_ranges[ncds_index:])

    elif transcript.strand == "-" and direction == "5":
        if cds_overlap:
            selected = list(cds_ranges[: cds_index + 1])
            s, e = selected[-1]
            selected[-1] = (pos, e)
        else:
            selected = list(cds_ranges[: ncds_index])

    else:  # strand == '-' and direction == '3'
        if cds_overlap:
            selected = list(cds_ranges[cds_index:])
            s, e = selected[0]
            selected[0] = (pos, e)
        else:
            selected = list(cds_ranges[ncds_index:])

    cut_length = 0
    for s, e in selected:
        if s is None or e is None or e < s:
            continue
        cut_length += (e - s + 1)

    if cut_length <= 0:
        return None

    return {"transcript": transcript, "part": direction, "cut_length": cut_length}


def cdscollection_nt_sequence(cc: dict) -> str:
    transcript = cc["transcript"]
    cut_length = cc["cut_length"]
    part = cc["part"]

    coding = transcript.coding_sequence
    comp_length = len(coding)
    if comp_length == 0 or cut_length <= 0:
        return ""

    if part == "5":
        start = 0
        end = min(comp_length - 1, start + cut_length - 1)
    else:
        end = comp_length - 1
        start = end - cut_length + 1 - 3
        if start < 0:
            start = 0

    if start > end:
        return ""
    return coding[start : end + 1]



def sv_transcript_build_dict(sv: dict, cc1: dict, cc2: dict):
    if not cc1 or not cc2:
        return None

    # Normalize insertion to a plain DNA string 
    ins = normalize_insseq(sv.get("insertion", ""))

    # insertion orientation depends on strand of cc1 transcript
    if cc1["transcript"].strand == "+":
        nt_ins = ins
    else:
        nt_ins = str(Seq(ins).reverse_complement())

    cc1_nt = cdscollection_nt_sequence(cc1)
    cc2_nt = cdscollection_nt_sequence(cc2)

    if cc1["part"] == "5":
        nt_cds = cc1_nt + nt_ins + cc2_nt
        nt_3utr = cc2["transcript"].three_prime_utr_sequence
    else:
        nt_cds = cc2_nt + nt_ins + cc1_nt
        nt_3utr = cc1["transcript"].three_prime_utr_sequence

    return {
        "sv": sv,
        "cc_1": cc1,
        "cc_2": cc2,
        "nt_sequence_ins": nt_ins,
        "nt_sequence_cds": nt_cds,
        "nt_sequence_3utr": nt_3utr,
    }



def sv_to_svtranscript(sv: dict, ensembl):
    transcript_1 = get_transcript(sv["chrom1"], sv["pos1"], ensembl, complete=True)
    transcript_2 = get_transcript(sv["chrom2"], sv["pos2"], ensembl, complete=True)
    if not transcript_1 or not transcript_2:
        return None

    cc1 = None
    cc2 = None
    p = sv["pattern"]

    if transcript_1.strand == "+" and transcript_2.strand == "+":
        if p == 1:
            cc1 = truncate_cds_simple(transcript_1, "5", sv["pos1"])
            cc2 = truncate_cds_simple(transcript_2, "3", sv["pos2"])
        elif p == 3:
            cc1 = truncate_cds_simple(transcript_1, "3", sv["pos1"])
            cc2 = truncate_cds_simple(transcript_2, "5", sv["pos2"])

    elif transcript_1.strand == "+" and transcript_2.strand == "-":
        if p == 2:
            cc1 = truncate_cds_simple(transcript_1, "5", sv["pos1"])
            cc2 = truncate_cds_simple(transcript_2, "3", sv["pos2"])
        elif p == 4:
            cc1 = truncate_cds_simple(transcript_1, "3", sv["pos1"])
            cc2 = truncate_cds_simple(transcript_2, "5", sv["pos2"])

    elif transcript_1.strand == "-" and transcript_2.strand == "+":
        if p == 2:
            cc1 = truncate_cds_simple(transcript_1, "3", sv["pos1"])
            cc2 = truncate_cds_simple(transcript_2, "5", sv["pos2"])
        elif p == 4:
            cc1 = truncate_cds_simple(transcript_1, "5", sv["pos1"])
            cc2 = truncate_cds_simple(transcript_2, "3", sv["pos2"])

    else:  # (-,-)
        if p == 1:
            cc1 = truncate_cds_simple(transcript_1, "3", sv["pos1"])
            cc2 = truncate_cds_simple(transcript_2, "5", sv["pos2"])
        elif p == 3:
            cc1 = truncate_cds_simple(transcript_1, "5", sv["pos1"])
            cc2 = truncate_cds_simple(transcript_2, "3", sv["pos2"])

    return sv_transcript_build_dict(sv, cc1, cc2)

def _five_three_cc(mut_transcript: dict):
    """Return (five_cc, three_cc) regardless of cc_1/cc_2 order."""
    cc1 = mut_transcript["cc_1"]
    cc2 = mut_transcript["cc_2"]
    if cc1["part"] == "5":
        return cc1, cc2
    return cc2, cc1


def is_inframe_junction(mut_transcript: dict) -> bool:
    """
    Canonical in-frame vs frameshift for a mut_transcript junction based on codon phase.
    Uses the actual lengths of the CDS fragments you assembled.
    """
    five_cc, three_cc = _five_three_cc(mut_transcript)

    five_nt = cdscollection_nt_sequence(five_cc)
    ins_nt = mut_transcript.get("nt_sequence_ins", "") or ""

    # phase after 5' + insertion
    phase_after_five = (len(five_nt) + len(ins_nt)) % 3

    # estimate the phase at the start of the 3' fragment in CDS coordinates
    three_tr = three_cc["transcript"]
    cds_len = len(three_tr.coding_sequence or "")
    three_cut = int(three_cc.get("cut_length", 0) or 0)

    # mirror cdscollection_nt_sequence(part="3") which uses: start = end - cut_length + 1 - 3
    three_start = max(0, cds_len - three_cut - 3)
    phase_at_three_start = three_start % 3

    return phase_after_five == phase_at_three_start


def frameshift_status(mut_transcript: dict) -> str:
    """
    Start-loss / In-frame / Frameshift.
    """
    if not mut_transcript or not mut_transcript.get("cc_1") or not mut_transcript.get("cc_2"):
        return "NA"

    five_cc, _three_cc = _five_three_cc(mut_transcript)
    five_nt = cdscollection_nt_sequence(five_cc)

    # 5' CDS fragment not starting at ATG
    if five_nt and not five_nt.startswith("ATG"):
        return "Start-loss"

    # phase-based junction check
    try:
        return "In-frame" if is_inframe_junction(mut_transcript) else "Frameshift"
    except Exception:
        return "NA"


# Sequence
# peptides + WT subtraction

def trim_to_3x(nt_sequence: str) -> str:
    r = len(nt_sequence) % 3
    return nt_sequence[:-r] if r else nt_sequence


def nt_to_aa(nt_sequence: str) -> str:
    dna_seq = Seq(trim_to_3x(nt_sequence))
    mrna_seq = dna_seq.transcribe()
    aa_seq = mrna_seq.translate(to_stop=True)
    return str(aa_seq)


def cut_sequence(aa_sequence: str, k: int):
    if not aa_sequence or len(aa_sequence) < k:
        return []
    return [aa_sequence[i:i+k] for i in range(len(aa_sequence) - k + 1)]


def mut_specific_peptides(mut_aa: str, wt_aas, k: int):
    """
    novel = mut_kmers - wt_kmers
    shared = mut_kmers & wt_kmers
    """
    mut_kmers = set(cut_sequence(mut_aa, k))
    wt_kmers = set()
    for wt in wt_aas:
        if wt:
            wt_kmers |= set(cut_sequence(wt, k))

    novel = sorted(mut_kmers - wt_kmers)
    shared = sorted(mut_kmers & wt_kmers)

    stats = {
        "mut_kmers": len(mut_kmers),
        "wt_kmers": len(wt_kmers),
        "novel_kmers": len(novel),
        "shared_kmers": len(shared),
    }
    return novel, shared, stats

def sv_gene_name(gene1, gene2, sep="_"):
    # normalize None 
    gene1 = (gene1 or "").strip()
    gene2 = (gene2 or "").strip()

    # if both empty
    if not gene1 and not gene2:
        return ""

    # if one empty
    if not gene1:
        return gene2
    if not gene2:
        return gene1

    # both non-empty
    if gene1 == gene2:
        return gene1
    return f"{gene1}{sep}{gene2}"



# Output

def write_tsv(results, out_path):
    fieldnames = [
        "ID",
        "sv_sorted_coord",
        "chrom1",
        "pos1",
        "chrom2",
        "pos2",
        "gene1", 
        "gene2",
        "sv_gene",
        "svlen",
        "pattern", 
        "svtype",
        "svtype_original",
        "somatic_score",
        "frameshift_status",
        "sv_gene_fusion_boolean",
        "mut_nt_sequence",
        "mut_aa_sequence",
        "wt_t1_aa_sequence",
        "wt_t2_aa_sequence",
        "k",
        "mut_kmers", 
        "wt_kmers", 
        "novel_kmers",
        "novel_kmers_sequence",
    ]

               
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()

        for r in results:
            row = dict(r)
            row["sv_sorted_coord"] = "|".join(row["sv_sorted_coord"]) if isinstance(row["sv_sorted_coord"], (list, tuple)) else row["sv_sorted_coord"]
            if isinstance(row.get("novel_kmers_sequence"), list):
                row["novel_kmers_sequence"] = ";".join(row["novel_kmers_sequence"])
            w.writerow({k: row.get(k, "") for k in fieldnames})



# Main pipeline


def run_summary(
    vcf_file: str,
    k: int = 9,
):
    #  standardized SV dicts 
    svs = vcf_load(vcf_file)
    # print(f'svs is {svs}')

    # dedup
    svs = remove_duplicate(svs)

    for sv in svs:
        if not isinstance(sv.get("insertion", ""), str):
            print("BAD insertion type:", type(sv["insertion"]), sv)
            break


    # load ensembl
    ensembl = ensembl_load("custom")

    # mut_transcript dict
    sv_transcripts = []
    for sv in svs:
        mut_transcript = sv_to_svtranscript(sv, ensembl)
        if mut_transcript:
            sv_transcripts.append(mut_transcript)

    # compute nt_sequence / aa_sequence / novel_kmers_sequence as VARIABLES in a loop
    results = []
    for i, mut_transcript in enumerate(sv_transcripts, start=1):
        nt_sequence = mut_transcript["nt_sequence_cds"] + mut_transcript["nt_sequence_3utr"]
        aa_sequence = nt_to_aa(nt_sequence)

        wt1 = getattr(mut_transcript["cc_1"]["transcript"], "protein_sequence", "") or ""
        wt2 = getattr(mut_transcript["cc_2"]["transcript"], "protein_sequence", "") or ""

        novel_kmers_sequence, shared_kmers, stats = mut_specific_peptides(
            mut_aa=aa_sequence,
            wt_aas=[wt1, wt2],
            k=k,
        )

        gene1 = get_gene_name(mut_transcript["sv"]["chrom1"], mut_transcript["sv"]["pos1"], ensembl)
        gene2 = get_gene_name(mut_transcript["sv"]["chrom2"], mut_transcript["sv"]["pos2"], ensembl)

        sv_gene = sv_gene_name(gene1, gene2)
        # # sv_gene rules:
        # # - if equal: gene1
        # # - else: gene1_gene2
        # if gene1 == gene2:
        #     sv_gene = gene1
        # else:
        #     sv_gene = f"{gene1}_{gene2}"

        sv_gene_fusion_boolean = bool(gene1) and bool(gene2) and (gene1 != gene2)
        fs = frameshift_status(mut_transcript)

        results.append(
            {
                "ID": f"sv{i}",
                "sv_sorted_coord": mut_transcript["sv"]["sorted_coord"],
                "chrom1": mut_transcript["sv"]["chrom1"],
                "pos1": mut_transcript["sv"]["pos1"],
                "chrom2": mut_transcript["sv"]["chrom2"],
                "pos2": mut_transcript["sv"]["pos2"],
                "gene1": gene1,
                "gene2": gene2,
                "sv_gene": sv_gene,
                "svlen": mut_transcript["sv"]["SVLEN"],
                "pattern": mut_transcript["sv"]["pattern"],
                "svtype": mut_transcript["sv"]["svtype"],
                "svtype_original": mut_transcript["sv"]["SVTYPE"],
                "somatic_score": mut_transcript["sv"]["SOMATICSCORE"],
                "frameshift_status": fs,
                "sv_gene_fusion_boolean": sv_gene_fusion_boolean,
                "mut_nt_sequence": nt_sequence,
                "mut_aa_sequence": aa_sequence,
                "wt_t1_aa_sequence": wt1,
                "wt_t2_aa_sequence": wt2,
                "k": k,
                "mut_kmers": stats["mut_kmers"],
                "wt_kmers": stats["wt_kmers"],
                "novel_kmers": stats["novel_kmers"],
                "novel_kmers_sequence": novel_kmers_sequence,
            }
        )

    return results



def main():
    ap = argparse.ArgumentParser(description="NeoSV pre-NetMHC pipeline (single file, dict-based)")
    ap.add_argument("--vcf", nargs="?", const=None, default=None, help="SV VCF file (.vcf or .vcf.gz)")
    ap.add_argument("--max_print", type=int, default=10, help="print first N mut_transcripts summary")
    ap.add_argument("--out", default="sv_peptides_results.tsv", help="output TSV file")
    ap.add_argument("--k", type=int, default=9, help="k-mer length for peptide generation (default 9)")
    args = ap.parse_args()

    # Treat None/"" as "no input"
    if not args.vcf:
        write_tsv([], args.out)
        print(f"No --vcf provided. Wrote header only to {args.out}")
        return

    # also treat non-existent file as "no input"
    if not os.path.exists(args.vcf):
        write_header_only([], args.out)
        print(f"VCF not found: {args.vcf}. Wrote header only to {args.out}")
        return

    results = run_summary(
        vcf_file=args.vcf,
        k=args.k,
    )

    write_tsv(results, args.out)

    print(f"\nDone. Produced {len(results)} mut_transcript results and wrote to {args.out}\n Printed the first {args.max_print} summary\n")
    for r in results[: args.max_print]:
        print(
            r["sv_sorted_coord"],
            "svtype=", r["svtype"],
            "AA_len=", len(r["mut_aa_sequence"]),
            "mut_kmers=", r["mut_kmers"],
            "novel_kmers=", r["novel_kmers"],
            "frameshift=", r["frameshift_status"],
        )


if __name__ == "__main__":
    main()

#example: python flatten_neosv_allSVTYPE.py --vcf SV_somatic_raw.manta.vcf.gz --out test_results_allSV.csv --k 9
