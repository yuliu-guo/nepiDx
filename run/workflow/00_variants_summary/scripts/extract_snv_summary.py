#!/usr/bin/env python3

#TODO: expand consequences list
#full list of valid_consequences = {"missense_variant", "frameshift_variant", "inframe_insertion", "inframe_deletion", 
#"stop_gained", "stop_lost", "start_lost", "protein_altering_variant", "splice_acceptor_variant",
#"splice_donor_variant", "splice_region_variant"}
# VEP sometimes uses multiple consequences joined with '&'
# cons_set = set((Consequence or "").split("&"))
# if cons_set.isdisjoint(valid_consequences):
#    continue



import os
import pysam
import gzip
import pandas as pd
import warnings
import vcfpy
import pileup_analysis
import extract_info_HGVSp
import generate_mt_fasta
import process_transcript_expression
import generate_frameshift_fasta
import fetch_external_db
import argparse
import datetime
import STARcounts_tpm
import time 
from collections import defaultdict


def parse_arguments():
    parser = argparse.ArgumentParser(description="Process VEP annotated VCF, BAM (optional), and transcript (optional) with RNA read base quality threshold (optional).")

    # Adding arguments
    parser.add_argument(
        "-v", "--vcf_path",
        type=str,
        required=True,
        help="Path to the VCF file."
    )
    
    parser.add_argument(
        "--reference_fasta",
        type=str,
        required=False,
        help="Path to the reference fasta file. Possible to read from config file"
    )


    parser.add_argument(
        "-b", "--bam_path",
        nargs="?",
        const=None,
        type=str,
        required=False,
        default=None,
        help="Path to the RNA BAM/CRAM file."
    )

    parser.add_argument(
        "-t", "--transcript_path",
        nargs="?",
        const=None,
        type=str,
        required=False,
        default=None,
        help="Path to the transcript file using RefSeq (*.refseq.quant.sf)."
    )

    parser.add_argument(
        "-s", "--STARcounts_path",
        nargs="?",
        const=None,
        type=str,
        required=False,
        default=None,
        help="Path to STAR counts file."
    )

    parser.add_argument(
        "-q", "--RNA_read_base_quality_threshold",
        type=int,
        default=30,
        help="RNA read base quality threshold (default: 30)."
    )
    
    parser.add_argument(
        "-p", "--patient_id",
        type=str,
        required=True,
        help="Patient ID for output directory and file naming."
    )
    parser.add_argument(
        "-o", "--output_path",
        type=str,
        required=True,
        help="Path to output path."
    )

    

    return parser.parse_args()
    

def read_vcf_to_dataframe(
    vcf_path, 
    reference_fasta,
    bam_path=None, 
    transcript_path=None, 
    STARcounts_path=None, 
    RNA_read_base_quality_threshold=30,
    profile=False):
    """
    main funtion to read vcf and bam files into a combined dataframe.

     Parameters:
    -----------
    vcf_path: str
        Path to VCF file.
    bam_path : str
        Path to corresponding BAM file.
    transcript_path: str
        Path to the transcript file.
    RNA_read_base_quality_threshold: int, optional
        Quality threshold for selecting supporting RNA read base.

    Returns:
    --------
    pd.DataFrame
        A dataframe containing variants summary.
        
    """
    # Optional profiler
    timings = defaultdict(float)
    def _tstart():
        return time.perf_counter()
        
    def _tadd(key, t0):
        timings[key] += (time.perf_counter() -t0)

    # Check up if file path is in correct form / exists
    def _is_ok_path(p):
        return bool(p and str(p).strip() and os.path.exists(p))
    
    bam_ok = _is_ok_path(bam_path)
    transcript_ok = _is_ok_path(transcript_path)
    star_ok = _is_ok_path(STARcounts_path)
    
    # One-time load of transcript expression lookup 
    transcript_lookup = None
    if transcript_ok:
      t0 = _tstart()
      transcript_lookup = process_transcript_expression.load_quant_lookup(transcript_path)
      _tadd("transcrip_tpm_load", t0)   
    
    # One-time load of STAR gene count expression lookup 
    star_tpm_lookup = None
    if star_ok:
      t0 = _tstart()
      star_tpm_lookup = STARcounts_tpm.load_star_tpm_lookup(STARcounts_path)
      _tadd("star_tpm_load", t0)  
    
    # One-time load for bam file
    bam_handle = None
    if bam_ok:
      t0 = _tstart()
      bam_handle = pileup_analysis.open_alignment(bam_path, reference_fasta)
      _tadd("bam_load", t0)  
    
    
    # Caches
    #protein_cache = {}
    transcript_cache = {}
    star_cache = {}
    pileup_cache = {}


    # Initiate a dataframe
    #df = pd.DataFrame()
    t0 = _tstart()
    #read vcf file
    reader = vcfpy.Reader.from_path(vcf_path)
    _tadd("read_vcf", t0)

    t0 = _tstart()
    #read header related to CSQ
    csq_info = reader.header.get_info_field_info("CSQ")
    # Extract the format
    csq_format = None
    if csq_info:
        description = csq_info.description
        # CSQ format is part of the description, after "Format: "
        if "Format:" in description:
            csq_format = description.split("Format: ")[1]
    csq_fields = csq_format.split('|')

    target_fields = ["Consequence", "IMPACT", "SYMBOL", "Feature_type", "Feature", "BIOTYPE", "ENSP", 
    "Protein_position", "Amino_acids", "HGVSp", "VARIANT_CLASS", "CDS_position", "Codons"]
    field_indices = {field: csq_fields.index(field) for field in target_fields if field in csq_fields}
    _tadd("csq_header_parse", t0)

    valid_consequences = {
        "frameshift_variant",
         "inframe_insertion",
          "inframe_deletion",
           "missense_variant"} 

    # Counter for unique IDs
    rows = []
    unique_id_counter = 1

    # Extract info from vcf
    for record in reader:
        #if not record.is_snv(): # only include snv (might be unnecessary)
        #    continue
        if 'PASS'not in record.FILTER: # only include passed variants
            continue
        #print(f'record is {record}')

        # Basic info
        FILTER = ",".join(str(i) for i in record.FILTER)
        chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alt = ",".join(alt.value for alt in record.ALT)
        QSS = record.INFO.get("QSS", None) #evidence that ALT frequency differs between tumor and normal
        SomaticEVS = record.INFO.get("SomaticEVS", None) #Somatic Empirical Variant Score (higher generally means more confident)
        
        #gnomADg_AF = record.INFO.get("gnomADg_AF", None)
        #sgt = record.INFO['SGT'] #SGT: Most likely somatic genotype excluding normal noise states
        
        ##### BAM pileup (fast path: one pass, no intermediate lists)
        RNA_VAF = RNA_ALT_N = RNA_TOTAL_N = None
        if bam_handle is not None:
            cache_key = (chrom, pos, ref, alt, RNA_read_base_quality_threshold)
            cached = pileup_cache.get(cache_key)
            if cached is None:
                t0 = _tstart()
                vaf, alt_n, total_n = pileup_analysis.pileup_counts_and_vaf(
                    bam_handle,
                    chrom,
                    pos,
                    ref,
                    alt,
                    quality_threshold=RNA_read_base_quality_threshold,
                )
                _tadd("bam_pileup_total", t0)
        
                RNA_VAF = round(vaf, 3)
                RNA_ALT_N = alt_n
                RNA_TOTAL_N = total_n
                pileup_cache[cache_key] = (RNA_VAF, RNA_ALT_N, RNA_TOTAL_N)
            else:
                RNA_VAF, RNA_ALT_N, RNA_TOTAL_N = cached
        
        # CSQ info
        csq_value = record.INFO.get("CSQ",[])
        if not csq_value:
            continue
        # Extract the values for each CSQ entry
        for entry in csq_value:
            fields = entry.split('|')
            entry_data = {}
            for field, index in field_indices.items():
                if len(fields) > index:
                    entry_data[field] = fields[index]
                else:
                    entry_data[field] = None  # Handle missing fields 
        
            Consequence = entry_data.get("Consequence")
            if Consequence not in valid_consequences:
                continue
            
            if entry_data.get("HGVSp") is None:
                continue
            ##### Protein mutation/peptides
            # Protein mutation/peptides
            
            if Consequence != "frameshift_variant":
                ensp = entry_data.get("ENSP")
                if not ensp:
                    continue
                t0 = _tstart()
                wt_seq = fetch_external_db.get_protein_seq(ensp)
                _tadd("protein_fetch_total", t0)
                if wt_seq is None:
                    # No sequence available for this ENSP (even after base-id normalization)
                    continue
                t0 = _tstart()
                flanked_mut = generate_mt_fasta.generate_mutated_protein_seq(
                    entry_data.get("HGVSp"), Consequence, wt_seq, flanking_sequence_length=8
                    )
                flanked_mut_52 = generate_mt_fasta.generate_mutated_protein_seq(
                    entry_data.get("HGVSp"), Consequence, wt_seq, flanking_sequence_length=52
                    )
                nmer = generate_mt_fasta.generate_nmer_peptides(
                    flanked_mut, nmer_length=9
                    )
                _tadd("mutate_and_nmer_total", t0)
            else:
                t0 = _tstart()
                flanked_mut = generate_frameshift_fasta.get_fs_mutated_protein(
                    Consequence,
                    entry_data.get("Feature"),
                    entry_data.get("VARIANT_CLASS"),
                    entry_data.get("CDS_position"),
                    entry_data.get("Codons"),
                    flanking_sequence_length=8)
                flanked_mut_52 = generate_frameshift_fasta.get_fs_mutated_protein(
                    Consequence,
                    entry_data.get("Feature"),
                    entry_data.get("VARIANT_CLASS"),
                    entry_data.get("CDS_position"),
                    entry_data.get("Codons"),
                    flanking_sequence_length=52)
                nmer = generate_mt_fasta.generate_nmer_peptides(
                    flanked_mut, nmer_length=9
                    )
                _tadd("frameshift_and_nmer_total", t0)
            

            ##### Transcript expression
            transcript_TPM = transcript_numreads = None
            
            if transcript_lookup is not None:
                feature = entry_data.get("Feature")
                if feature:
                    t0 = _tstart()
                    vals = transcript_lookup.get(feature)  # (TPM, NumReads) or None
                    _tadd("transcript_lookup_total", t0)
            
                    if vals is not None:
                        tpm, nr = vals
                        transcript_TPM = round(tpm, 3) if tpm is not None else None
                        transcript_numreads = round(nr, 3) if nr is not None else None
                        
            ##### Gene expression
            gene_count_TPM = None
            if star_tpm_lookup is not None:
                t0 = _tstart()
                symbol = entry_data.get("SYMBOL")
                if symbol:
                    val = star_tpm_lookup.get(symbol)
                    gene_count_TPM = round(val, 3) if val is not None else None
                    _tadd("star_lookup_total", t0)

            unique_id = f"s{unique_id_counter}"
            unique_id_counter += 1

            rows.append({
                "ID": unique_id,
                "CHROM": chrom,
                "POS": pos,
                "REF": ref,
                "ALT": alt,
                "QSS": QSS,
                "SomaticEVS": SomaticEVS,
                "FILTER": FILTER,
                "Consequence": Consequence,
                "IMPACT": entry_data.get("IMPACT"),
                "SYMBOL": entry_data.get("SYMBOL"),
                "refseq": entry_data.get("Feature"),
                "ENSP": entry_data.get("ENSP"),
                "HGVSp": entry_data.get("HGVSp"),
                "Protein_position": entry_data.get("Protein_position"),
                "Amino_acids": entry_data.get("Amino_acids"),
                "transcript_TPM": transcript_TPM,
                "transcript_NumReads": transcript_numreads,
                "gene_count_TPM": gene_count_TPM,
                "flanked_mutated_protein_sequence": flanked_mut,
                "flanked_mutated_protein_sequence_52aa": flanked_mut_52,
                "nmer": nmer,
                "RNA_ALT_N": RNA_ALT_N,
                "RNA_TOTAL_N": RNA_TOTAL_N,
                "RNA_VAF": RNA_VAF,
            })

    t0 = _tstart()
    df = pd.DataFrame(rows)
    _tadd("dataframe_build", t0)
    
    #close opened files
    if bam_handle is not None:
      bam_handle.close()
    

    if profile:
        print("=== Timing breakdown (seconds) ===")
        for k, v in sorted(timings.items(), key=lambda kv: -kv[1]):
            print(f"{k:28s} {v:10.3f}")
        print("=== Inputs enabled ===")
        print(f"BAM enabled:        {bam_ok}")
        print(f"Transcript enabled: {transcript_ok}")
        print(f"STAR enabled:       {star_ok}")

    return df


    
if __name__ == "__main__":
    t0 = time.perf_counter()
    # Parse arguments
    args = parse_arguments()

    # Extract arguments
    vcf_file_path = args.vcf_path
    reference_fasta = args.reference_fasta
    bam_file_path = args.bam_path
    transcript_path = args.transcript_path
    STARcounts_path = args.STARcounts_path
    quality_threshold = args.RNA_read_base_quality_threshold
    patient_id = args.patient_id
    output_path = args.output_path
    
    if not os.path.exists(output_path):
        # Load the VCF data into a DataFrame


        df = read_vcf_to_dataframe(vcf_file_path, reference_fasta, bam_file_path, transcript_path, STARcounts_path, quality_threshold, profile=True)
        elapsed = time.perf_counter() - t0
        print(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} VCF data loaded successfully for patient {patient_id}.")

        # Save the DataFrame to an Excel file
        df.to_csv(output_path, sep='\t', index=False)
        elapsed = time.perf_counter() - t0
        print(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Data saved to {output_path}")
        print(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Total runtime: {elapsed:.2f} seconds")

    else:
        # Read the DataFrame from the existing Excel file
        df = pd.read_csv(output_path, sep='\t')
        elapsed = time.perf_counter() - t0
        print(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Data loaded from existing file {output_path}")
        print(f"{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Total runtime: {elapsed:.2f} seconds")

