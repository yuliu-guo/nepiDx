import re
from typing import Optional, Dict, Any

"""
This script is to support extracting information from HGVSp column in VEP annotated vcf file.

"""


def separate_letters_numbers(s):
    # Find all sequences of letters and numbers in the string
    letter_sequences = re.findall(r'[A-Za-z]+', s)
    number_sequences = re.findall(r'\d+', s)
    # e.g. Leu287
    if len(letter_sequences) == 1 and len(number_sequences) == 1:
        aa_start = letter_sequences[0]
        position_start = number_sequences[0]
        position_end = number_sequences[0]
        aa_end = letter_sequences[0]
        return {
            "aa_start":aa_start, 
            "aa_end": aa_end, 
            "position_start": position_start, 
            "position_end": position_end
            }

    # If there are multiple sequences of numbers and letters, process them
    # e.g. Ala496_Ala502
    # e.g. Leu20_Leu21
    if len(letter_sequences) == 2 and len(number_sequences) == 2:
        aa_start = letter_sequences[0]
        position_start = number_sequences[0]
        position_end = number_sequences[1]
        aa_end = letter_sequences[1]
        return {
            "aa_start":aa_start, 
            "aa_end": aa_end, 
            "position_start": position_start, 
            "position_end": position_end
            }
    # two amino acid, one postion #e.g. Ala124Val
    if len(letter_sequences) == 2 and len(number_sequences)== 1: 
        WT_aa = letter_sequences[0] if letter_sequences else ''
        MT_aa = letter_sequences[-1] if letter_sequences else ''
        position = ''.join(number_sequences)
        return {
            "WT_aa":WT_aa, 
            "MT_aa": MT_aa, 
            "position": position
            }
    else:
        print(f'more than 2 postions. Please specify')


def aa_3_to_1(aa_sequence):
    """Convert three-letter amino acid code to one-letter code for a concatenated sequence."""
    d = {
        'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K', 
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
        'TER': '*', 'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 
        'MET': 'M', 'XAA': 'X'
    }
    one_letter_code = ''
    for i in range(0, len(aa_sequence), 3):
        triplet = aa_sequence[i:i+3].upper()
        if triplet in d:
            one_letter_code += d[triplet]
        else:
            # Handle invalid triplet or add custom error handling
            one_letter_code += '?'
            print(f"Warning: Triplet '{triplet}' not found in dictionary.")
    return one_letter_code


def extract_parts(info):
    return separate_letters_numbers(info)

def process_missense_variant(info):
    '''
    with missense variant, the format is always WT_aa+position+MT_aa
    e.g. Ala124Val
    Ala is the wild type amino acid
    124 is the postion where the mutation happend
    Val is the mutated amino acid


    '''
    WT_aa_variable = None
    MT_aa_variable = None
    position_variable = None 
    
    parts = extract_parts(info)
    WT_aa_variable = aa_3_to_1(parts["WT_aa"])
    MT_aa_variable = aa_3_to_1(parts["MT_aa"])
    position_variable = int(parts["position"])
    return {
        "WT_aa": WT_aa_variable, 
        "MT_aa": MT_aa_variable, 
        "position": position_variable
    }


def process_inframe_deletion(info):
    """
    e.g. 
    Ala496_Ala502del
    Glu740del
    """
    aa_start_variable = None
    aa_end_variable = None
    position_start_variable = None
    position_end_variable = None

    variant_info = info.split("del")[0] # seperate variant info before 'del'
    parts = extract_parts(variant_info)
    aa_start_variable = aa_3_to_1(parts["aa_start"])
    aa_end_variable = aa_3_to_1(parts["aa_end"])
    position_start_variable = int(parts["position_start"])
    position_end_variable = int(parts["position_end"])
    return {
        "aa_start": aa_start_variable, 
        "aa_end": aa_end_variable, 
        "position_start": position_start_variable,
        "position_end": position_end_variable
    }

def process_inframe_insertion(info):
    #NP_001291362.1:p.Pro506_Cys507insSer
    #NP_001359162.1:p.Ala336dup
    aa_start_variable = None
    aa_end_variable = None
    position_start_variable = None
    position_end_variable = None
    ins_aa = None
    variant_info = None

    if "dup" in info:
        variant_info = info.split("dup")[0]
        ins_aa = ""
    elif "ins" in info:
        if 'del' in info:
            print(f'Inframe insersion {info} has irregular form - check VCF file. ')
            return None
        else:
            variant_info = info.split("ins")[0]
            ins_aa = aa_3_to_1(info.split("ins")[1])
# exception : NP_001180264.1:p.Phe66delinsLeuMetIle
# variant_info = Phe66del #Phe = F
# ins_aa = LMI

    parts = extract_parts(variant_info)
    aa_start_variable = aa_3_to_1(parts["aa_start"])
    aa_end_variable = aa_3_to_1(parts["aa_end"])
    position_start_variable = int(parts["position_start"])
    position_end_variable = int(parts["position_end"])
    return {
        "aa_start": aa_start_variable, 
        "aa_end": aa_end_variable, 
        "position_start": position_start_variable,
        "position_end": position_end_variable,
        "ins_aa": ins_aa
    }

def process_HGVSp(HGVSp: Optional[str], Consequence: str) -> Optional[Dict[str, Any]]:
    """
    HGVSp: str or None
    """
    if HGVSp is None or ":p." not in HGVSp:
        return None    
    else:
        p_info = HGVSp.split(":p.")[1]

    if Consequence == "missense_variant":
        return process_missense_variant(p_info)
    elif Consequence == "frameshift_variant":
        return process_frameshift_variant(p_info)
    elif Consequence == "inframe_deletion":
        return process_inframe_deletion(p_info)
    elif Consequence == "inframe_insertion":
        return process_inframe_insertion(p_info)
    else:
        return None

#print(process_HGVSp("NP_001180264.1:p.Phe66delinsLeuMetIle","inframe_insertion"))
