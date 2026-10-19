import os
import re
import extract_info_HGVSp

"""
This script is to generate mutated peptides from provided information from HGVSp column in VEP annotated vcf file

Current supporting variants are: missense_variant, inframe_deletion, inframe_insertion

Note: Frameshift variants are processed in a seperate script


"""

def get_protein_seq(ensp_id):
    fastaFile = "../ref/GRCh38_latest_protein.faa"
    '''Fetches the protein sequence for a given ENSP ID from a fasta file'''
    with open(fastaFile, 'r') as file:
        sequence_id = None
        sequence = []
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                if sequence_id == ensp_id:
                    return ''.join(sequence)
                sequence_id = line[1:].split()[0]
                sequence = []
            else:
                if sequence_id == ensp_id:
                    sequence.append(line)
        if sequence_id == ensp_id:
            return ''.join(sequence)
    return None


def generate_missense_variant_protein_seq(HGSVp, Consequence, wt_protein_sequence, flanking_sequence_length=8):
    # check if consequence matches with the function purpose
    # in this case, only missense_variant should be include
    if Consequence != "missense_variant":
        print(f"Warning: funtion generate_missense_variant_protein_seq doesn't support consequence:{Consequence}")

    # process HGVSp information
    results = extract_info_HGVSp.process_HGVSp(HGVSp=HGSVp, Consequence=Consequence)
    if not results or not isinstance(results, dict):
        print(f"Failed to process HGVSp: {HGSVp} with Consequence: {Consequence}")
        flanked_mutated_protein_sequence = None
        flanked_wt_protein_sequence = None
    else:
        # Consequence based - missense_variant share the same format but different than other forms of variants
        WT_aa = results["WT_aa"] #get the info wildtype amino acid from HGVSp
        MT_aa = results["MT_aa"] #get the mutated amino acid
        position = results["position"] - 1 # get the amino acid position information from HGVSp and convert to 0-based index
        position_start = position_end = position  # missense_vaiant only have one mutation position provided
            
        
        # when wildtype sequence cannot be found based on provided ensp_id, sequence = None
        # check if amino acid is within the range of wildtype protein sequence
        if wt_protein_sequence != None and 0 <= position < len(wt_protein_sequence):
            #check if WT_aa matchs with sequence[pos_start]
            if wt_protein_sequence[position] != WT_aa:
                print(f'WARNING: {HGSVp}: {Consequence}: Original amino acid {wt_protein_sequence[position]} does not match expected {WT_aa} at position {position+1}')
            
            # convert protein sequence to a list and replace with mutated amino acid at given postion
            wt_protein_sequence = list(wt_protein_sequence)
            wt_protein_sequence[position] = MT_aa
            mt_protein_sequence = ''.join(wt_protein_sequence)
            
            # extract the mutated amino acid with flanking on both side 
            # defalut flanking_sequence_length is set to 8
            # In case the flanking start/end exceeded the wildtype length
            flanking_position_start = max(position_start - flanking_sequence_length , 0) 
            flanking_position_end = min(position_end + flanking_sequence_length, len(wt_protein_sequence))

            flanked_mutated_protein_sequence = mt_protein_sequence[flanking_position_start:flanking_position_end]
            
        else:
            print(f'WARNING:{HGSVp}: {Consequence}: Position {position} is out of range for Sequence {wt_protein_sequence}')
            flanked_mutated_protein_sequence = None
           

    return flanked_mutated_protein_sequence

def generate_inframe_deletion_variant_protein_seq(HGSVp, Consequence, wt_protein_sequence, flanking_sequence_length=8):

    if Consequence != "inframe_deletion":
        print(f"funtion generate_inframe_deletion_variant_protein_seq doesn't support consequence:{Consequence}")

    results = extract_info_HGVSp.process_HGVSp(HGVSp=HGSVp, Consequence=Consequence)
    if not results or not isinstance(results, dict):
        print(f"Failed to process HGVSp: {HGSVp} with Consequence: {Consequence}")
        print(f'more info: {results}')
        flanked_mutated_protein_sequence = None
        #flanked_wt_protein_sequence = None
    
    else: 
        aa_start = results["aa_start"]
        aa_end = results["aa_end"]
        position_start = results["position_start"] - 1
        position_end= results["position_end"] - 1

        if wt_protein_sequence != None:
            #check if aa_start/end match with sequence[pos_start/end]
            if wt_protein_sequence[position_start] != aa_start:
                print(f"Amino acid at position {position_start} is {wt_protein_sequence[position_start]} and does not match with {aa_start} ")
            if wt_protein_sequence[position_end] != aa_end:
                print(f"Amino acid at position {position_end} is {wt_protein_sequence[position_end]} and does not match with {aa_end} ")
            # replace mutated amino acid and get mutated protein sequence
            if position_start < position_end:
                mt_protein_sequence = wt_protein_sequence[:position_start]+ wt_protein_sequence[position_end:]
            if position_start == position_end:
                position_end = position_end + 1
                mt_protein_sequence = wt_protein_sequence[:position_start]+ wt_protein_sequence[position_end:]

            # get flanked mutated sequence
            flanking_position_start = max(position_start - flanking_sequence_length , 0) 
            flanking_position_end = min(position_end + flanking_sequence_length, len(wt_protein_sequence))

            flanked_mutated_protein_sequence = mt_protein_sequence[flanking_position_start:flanking_position_end]
            
        else:
            print(f'Position {position_start}-{position_end} is out of range for Sequence {wt_protein_sequence}, Consequence is {Consequence}')
            flanked_mutated_protein_sequence = None

    return flanked_mutated_protein_sequence


def generate_inframe_insertion_variant_protein_seq(HGSVp, Consequence, wt_protein_sequence, flanking_sequence_length=8):

    if Consequence != "inframe_insertion":
        print(f"funtion generate_inframe_insertion_variant_protein_seq doesn't support consequence:{Consequence}")
        flanked_mutated_protein_sequence = None
        flanked_wt_protein_sequence = None

    results = extract_info_HGVSp.process_HGVSp(HGVSp=HGSVp, Consequence=Consequence)
    if not results or not isinstance(results, dict):
        print(f"Failed to process HGVSp: {HGSVp} with Consequence: {Consequence}")
        flanked_mutated_protein_sequence = None
        #flanked_wt_protein_sequence = None

    else:
        aa_start = results["aa_start"]
        aa_end = results["aa_end"]
        position_start = results["position_start"] - 1
        position_end= results["position_end"] - 1
        ins_aa = results["ins_aa"]

        if wt_protein_sequence != None:
            #check if aa_start/end match with sequence[pos_start/end]
            if wt_protein_sequence[position_start] != aa_start:
                print(f"Amino acid at position {position_start} is {wt_protein_sequence[position_start]} and does not match with {aa_start} ")
            if wt_protein_sequence[position_end] != aa_end:
                print(f"Amino acid at position {position_end} is {wt_protein_sequence[position_end]} and does not match with {aa_end} ")
            
            if "dup" in HGSVp:
                if position_start < position_end:
                    ins_aa = wt_protein_sequence[position_start:position_end]
                    mt_protein_sequence = wt_protein_sequence[:position_start]+ ins_aa + wt_protein_sequence[position_end:]
                
                if position_start == position_end:
                    ins_aa = wt_protein_sequence[position_start]
                    mt_protein_sequence = wt_protein_sequence[:position_start]+ ins_aa + wt_protein_sequence[position_end:]

            if "ins" in HGSVp:
                mt_protein_sequence = wt_protein_sequence[:position_start+1]+ ins_aa + wt_protein_sequence[position_end:]

            # get flanked mutated sequence
            flanking_position_start = max(position_start - flanking_sequence_length , 0) 
            flanking_position_end = min(position_end + flanking_sequence_length, len(wt_protein_sequence))

            flanked_mutated_protein_sequence = mt_protein_sequence[flanking_position_start:flanking_position_end]
            
        else:
            print(f'Position {position_start}-{position_end} is out of range for Sequence {wt_protein_sequence}, Consequence is {Consequence}')
            flanked_mutated_protein_sequence = None

    return flanked_mutated_protein_sequence

def generate_mutated_protein_seq(HGSVp, Consequence, wt_protein_sequence, flanking_sequence_length=8):
    if Consequence == "missense_variant":
        return generate_missense_variant_protein_seq(HGSVp=HGSVp, 
        Consequence=Consequence, 
        wt_protein_sequence=wt_protein_sequence, 
        flanking_sequence_length=flanking_sequence_length)
    elif Consequence == "inframe_deletion":
        return generate_inframe_deletion_variant_protein_seq(HGSVp=HGSVp, 
        Consequence=Consequence, 
        wt_protein_sequence=wt_protein_sequence, 
        flanking_sequence_length=flanking_sequence_length)
    elif Consequence == "inframe_insertion":
        return generate_inframe_insertion_variant_protein_seq(HGSVp=HGSVp, 
        Consequence=Consequence, 
        wt_protein_sequence=wt_protein_sequence, 
        flanking_sequence_length=flanking_sequence_length)

    else:
        print(f'Consequence {Consequence} is not supported')
        return None


def generate_nmer_peptides(flanked_mutated_protein_sequence, nmer_length=9):
    """
    with provided flanked mutated peptide, do sliding window to generate ideal length of n_mer peptides

    return a list of nmer peptides (or an empty list)
    """
    nmers = []
    if flanked_mutated_protein_sequence != None: 
    #if position_start is not None and position_end is not None:
        for i in range(len(flanked_mutated_protein_sequence)):
            if i >= 0 and i + nmer_length <= len(flanked_mutated_protein_sequence):
                nmers.append(''.join(flanked_mutated_protein_sequence[i:i + nmer_length]))
    return nmers


