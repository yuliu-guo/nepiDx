# make this as an option function

import vcfpy
import csv
import sys
import os
from abc import ABCMeta
from collections import OrderedDict
import binascii
import re
import glob
from Bio import SeqIO
import logging
import pandas as pd
import get_sample_id
import argparse

def parse_arguments():
    parser = argparse.ArgumentParser(description="read path and return only the sample id")

    parser.add_argument(
        "--arriba_file",
        type=str,
        required=True,
        help="Full path to the Arriba output file."
    )
    parser.add_argument(
        "--starfusion_file",
        nargs="?",
        const=None,
        type=str,
        required=False,
        default=None,
        help="Full path to the STAR-Fusion output file (optional)."
    )

    parser.add_argument(
        "--fusionmap_file",
        nargs="?",
        const=None,
        type=str,
        required=False,
        default=None,
        help="Full path to the FusionMap output file (optional)."
    )
    parser.add_argument(
        "-o", "--output_file",
        type=str,
        required=True,
        help="output file path"
    )
    return parser.parse_args()


def construct_index(count, gene, transcript, variant_type, position):
    return '{}.{}.{}.{}.{}'.format(count, gene, transcript, variant_type, position)

def determine_fusion_sequence(full_sequence, separator):
    if separator not in full_sequence:
        sys.exit("Fusion position marker '*' not found in fusion sequence. If running with AGFusion results, please rerun AGFusion using the `--middlestar` option.")
    else:
        fusion_position = full_sequence.find(separator)
        sequence = full_sequence.replace(separator, '')
        return (fusion_position, sequence)

def parse_starfusion_file(starfusion_file):
    if not os.path.exists(starfusion_file):
        raise Exception("Starfusion file {} not found. Aborting.".format(starfusion_file))
    with open(starfusion_file, 'r') as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)

def parse_fusionmap_file(fusionmap_file):
    if not os.path.exists(fusionmap_file):
        raise Exception ("FusionMap file {} not found. Aborting.".format(fusionmap_file))
    with open(fusionmap_file, 'r') as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)

def find_matching_starfusion_entry_for_arriba(starfusion_entries, five_prime_chr, five_prime_start, three_prime_chr, three_prime_start):
    for starfusion_entry in starfusion_entries:
        (sf_five_prime_chr, five_prime_pos, five_prime_strand) = starfusion_entry['LeftBreakpoint'].split(':')
        (sf_three_prime_chr, three_prime_pos, three_prime_strand) = starfusion_entry['RightBreakpoint'].split(':')
        if (sf_five_prime_chr.replace("chr", "") == five_prime_chr.replace("chr", "")
            and sf_three_prime_chr.replace("chr", "") == three_prime_chr.replace("chr", "")):
            if five_prime_pos == five_prime_start and three_prime_pos == three_prime_start:
                return starfusion_entry

def find_matching_fusionmap_filtered_entry_for_arriba(fusionmap_entries, five_prime_chr, five_prime_start, three_prime_chr, three_prime_start):
    for fusionmap_entry in fusionmap_entries:
        sf_five_prime_chr = fusionmap_entry['Chromosome1']
        five_prime_pos = fusionmap_entry['Position1'] 
        sf_three_prime_chr = fusionmap_entry['Chromosome2']
        three_prime_pos = fusionmap_entry['Position2']
        if (sf_five_prime_chr.replace("chr", "") == five_prime_chr.replace("chr", "")
            and sf_three_prime_chr.replace("chr", "") == three_prime_chr.replace("chr", "")):
            if five_prime_pos == five_prime_start and three_prime_pos == three_prime_start:
                return fusionmap_entry

def parse_arriba_file(input_file, starfusion_entries, fusionmap_entries):
        if not os.path.exists(input_file):
            raise Exception("Input file {} doesn't exist. Please provide a valid Arriba result file path. Aborting.".format(self.input_file))
        # Counter for unique IDs
        unique_id_counter = 1
        output_rows = []
        count = 1
        with open(input_file, 'r') as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for record in reader:
                (five_prime_chr, five_prime_start) = record['breakpoint1'].split(':')
                (three_prime_chr, three_prime_start) = record['breakpoint2'].split(':')
                if record['peptide_sequence'] == '.' or record['peptide_sequence'] is None or record['peptide_sequence'] == "":
                    continue
                (fusion_position, fusion_amino_acid_sequence) = determine_fusion_sequence(record['peptide_sequence'], '|')
                gene_name = "{}-{}".format(record['#gene1'], record['gene2'])
                transcript_name = "{}-{}".format(record['transcript_id1'], record['transcript_id2'])
                if record['reading_frame'] == 'in-frame':
                    variant_type = 'inframe_fusion'
                elif record['reading_frame'] == 'out-of-frame':
                    variant_type = 'frameshift_fusion'
                else:
                    continue
                if starfusion_entries is not None:
                    starfusion_entry = find_matching_starfusion_entry_for_arriba(starfusion_entries, five_prime_chr, five_prime_start, three_prime_chr, three_prime_start)
                else:
                    starfusion_entry = None
                if fusionmap_entries is not None:
                    fusionmap_entry = find_matching_fusionmap_filtered_entry_for_arriba(fusionmap_entries, five_prime_chr, five_prime_start, three_prime_chr, three_prime_start)
                else:
                    fusionmap_entry = None
                
                 # Add unique ID for each row
                unique_id = f"g{unique_id_counter}"
                unique_id_counter += 1

                output_row = {
                    "ID"                         : unique_id,
                    "CHROM"                      : "{} / {}".format(five_prime_chr, three_prime_chr),
                    "POS"                        : "{} / {}".format(five_prime_start, three_prime_start),
                    "REF"                        : 'fusion',
                    "ALT"                        : 'fusion',
                    "Protein_position"           : fusion_position,
                    "gene_name"                  : gene_name,
                    "variant_type"               : variant_type,
                    "fusion_amino_acid_sequence" : fusion_amino_acid_sequence.replace("*", "").upper(),
                    "transcript_name"            : transcript_name,
                    "index"                      : construct_index(count, gene_name, transcript_name, variant_type, fusion_position),
                    "fusion_read_support"        : int(record['split_reads1']) + int(record['split_reads2']) + int(record['discordant_mates']),
                }
                if starfusion_entry is not None:
                    output_row['fusion_expression']   = starfusion_entry['FFPM']
                    output_row['StarFusion']   = 'yes'
                else:
                    output_row['fusion_expression'] = 'NA'
                    output_row['StarFusion']   = 'no'

                if fusionmap_entry is not None:
                    output_row['fusionMap'] = 'yes'
                else:
                    output_row['fusionMap'] = 'no'
                


                output_rows.append(output_row)
                count += 1
        return output_rows


def main(input_file, output_file, starfusion_file=None, fusionmap_file=None):
    # Parse StarFusion file if provided
    starfusion_entries = None
    if starfusion_file:
        starfusion_entries = parse_starfusion_file(starfusion_file)
    
    fusionmap_entries = None
    if fusionmap_file:
        fusionmap_entries = parse_fusionmap_file(fusionmap_file)


    # Process Arriba file or directory
    if os.path.isfile(input_file):
        output_rows = parse_arriba_file(input_file, starfusion_entries, fusionmap_entries)
    else:
        raise Exception(f"Input file {input_file} doesn't exist. Aborting.")
    
    #convert output rows to a dataframe
    df = pd.DataFrame(output_rows)

    # write empty file with headers
    if df.empty: 
        df.to_csv(output_file, sep="\t", index=False)
        return

    # Only filter if at least one external caller file was provided
    if starfusion_file or fusionmap_file:
        # if column missing for some reason, create it as empty
        if "StarFusion" not in df.columns:
            df["StarFusion"] = ""
        if "fusionMap" not in df.columns:
            df["fusionMap"] = ""

        df = df[(df["StarFusion"] == "yes") | (df["fusionMap"] == "yes")]


    # Write output rows to the file
    #df.to_csv(output_file, index=False,sep='\t')
    df.to_csv(output_file, sep = '\t', index=False)


def output_headers():
        return[
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "Protein_position",
        "gene_name",
        "variant_type",
        "fusion_amino_acid_sequence",
        "transcript_name",
        "index",
        "fusion_read_support",
        "fusion_expression",
        "StarFusion",
        "fusionMap",

        ]


    
if __name__ == "__main__":
    # Parse arguments
    args = parse_arguments()

    # Extract arguments
    #file_path = args.path
    arriba_file = args.arriba_file
    starfusion_file = args.starfusion_file
    fusionmap_file = args.fusionmap_file
    #prefix = args.prefix
    output_file = args.output_file
    

if starfusion_file and fusionmap_file:
    main(arriba_file, output_file, starfusion_file=starfusion_file, fusionmap_file=fusionmap_file)
elif starfusion_file and not fusionmap_file:
    main(arriba_file, output_file, starfusion_file=starfusion_file)
elif (not starfusion_file) and fusionmap_file:
    main(arriba_file, output_file, fusionmap_file=fusionmap_file)
else:
    main(arriba_file, output_file)
