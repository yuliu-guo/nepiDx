import os
import re
from Bio.Seq import Seq
from Bio.Data import CodonTable
import extract_info_HGVSp
import fetch_external_db
import math


### a combined function to get CCDS sequence from RefSeq transcript ID
def get_CCDS_seq_from_RefSeq(RefSeq_transcript_id):
    '''Fetches the CCDS sequence for a given RefSeq transcript ID'''
    CCDS_id = fetch_external_db.get_CCDS_ID(RefSeq_transcript_id) # get CCDS ID from RefSeq transcript ID
    if CCDS_id is not None:
        return fetch_external_db.get_CCDS_seq(CCDS_id) # get CCDS sequence from CCDS ID
    return None

def translate_CCDS_with_start_stop(CCDS_seq):
    #standard_table = CodonTable.unambiguous_dna_by_name["Standard"]
    
    # Find the first occurrence of a start codon ("ATG" in DNA)
    start_pos = CCDS_seq.find("ATG")
    
    if start_pos == -1:
        return "No start codon found"
    else:
        # Extract the sequence starting from the start codon
        seq_from_start = CCDS_seq[start_pos:]
        #seq_from_start = Seq(CCDS_seq[start_pos:])
        # Add trailing 'N' to make the length a multiple of three 
        # (disabled due to unwanted generation of ambiguous amino acid "X")
        #while len(seq_from_start) % 3 != 0:
        #  seq_from_start += 'N'

        #Trim trailing bases 
        trim_len = len(seq_from_start) - (len(seq_from_start) % 3)
        seq_from_start = seq_from_start[:trim_len]
    
        # Translate the sequence, stopping at the first stop codon
        protein_seq = Seq(seq_from_start).translate(table="Standard", to_stop=True)
    
    return protein_seq

def get_fs_mutated_protein(Consequence, Feature, VARIANT_CLASS, CDS_position, Codons, flanking_sequence_length=8):
    '''Fetches the CDS sequence for a given Feature, VARIANT_CLASS, CDS_position, Codons and returns mutated protein information.'''
    if Consequence != "frameshift_variant":
        print(f"Variant consequence {Consequence} is not supported by this function get_fs_mutated_protein")
    
    # first get CCDS sequence from nucleotide id, which is recorded in Feature.
    CCDS_seq = get_CCDS_seq_from_RefSeq(Feature)

    if CCDS_seq != None and CDS_position != None and VARIANT_CLASS != None: 
        # second, get the mutation position from CDS_postion, including the start and end position and converted to 0 based. 
        if CDS_position != "":
            if '-' in CDS_position:
                start_pos = int(CDS_position.split('-')[0]) - 1
                end_pos = int(CDS_position.split('-')[1]) - 1
            else: 
                start_pos = int(CDS_position) - 1
                end_pos = start_pos
        else:
            print("Invalid CDS_position format")
        
        approx_start_pos_on_aa_level = math.floor(start_pos/3) # NOTE: floor works for 0-based, ceil works for 1-based

        if VARIANT_CLASS == "insertion":
            #if len(ref) < len(alt): #(two ways of identifying insertion/deletion)
            # If the variant is an insertion, the CDS position refers to the base before the insertion
            # Codon format: -(ref)/GGCTGC(alt) for insertion
            # Codon format: GGCTGC(ref)/-(alt) for deletion
            if Codons != None:
                codon = Codons.split('/')[1] #take the inserted codon
                extended_aa_len = round(len(codon)/3) # estimate the length of mutated amino acid 
        
                #get the full CCDS sequence, which used for further translation
                mutated_CCDS_seq = CCDS_seq[:start_pos] + codon + CCDS_seq[start_pos:] 
            else:
                print(f"Codons doesn't exist for {Feature} and {VARIANT_CLASS}")
        
        elif VARIANT_CLASS == "deletion":
        #elif len(ref) > len(alt): #(two ways of identifying insertion/deletion)
            # If the variant is a deletion, the CDS position refers to the first base that's deleted
            mutated_CCDS_seq = CCDS_seq[:start_pos] + CCDS_seq[end_pos + 1:]
            extended_aa_len = 0
        else:
            raise ValueError("Unsupported VARIANT_CLASS")

        # finally get the mutated translated protein sequence
        mt_protein_sequence = translate_CCDS_with_start_stop(mutated_CCDS_seq)

        # adding flanking 
        # get flanked mutated sequence
        flanking_position_start = max(approx_start_pos_on_aa_level - flanking_sequence_length , 0) 
        #until the end of the mutated protein because noramlly from the starting position, all the downstream Aa are also changed due to frameshift
        flanked_mutated_protein_sequence = mt_protein_sequence[flanking_position_start:]
    else:
        return None
    return flanked_mutated_protein_sequence

# Test the function
#Feature = "NM_001127208.3"
#VARIANT_CLASS = "insertion"
#CDS_position = "3813"
#CDS_position = "3813-3814"
#Codons = "-/G"
#Codons = "-/GCCTGTCAGGNN"
#HGSVp = "NP_001120680.1:p.Asp1277SerfsTer26"
#Consequence="frameshift_variant"

#print(get_fs_mutated_protein(Consequence, Feature, VARIANT_CLASS, CDS_position, Codons, 9))
