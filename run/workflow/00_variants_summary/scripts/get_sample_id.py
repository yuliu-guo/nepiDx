import re
from pathlib import Path
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(description="read path and return only the sample id")

    # Adding arguments
    parser.add_argument(
        "-p", "--path",
        type=str,
        required=True,
        help="Path to VCF/CRAM/BAM or other file."
    )
    parser.add_argument(
        "-pf", "--prefix",
        type=str,
        required=True,
        help="Prefix of the sample id, it can be Path to replace the original path."
    )
    
    parser.add_argument(
        "-sf", "--suffix",
        type=str,
        required=True,
        help="suffix of the sample id, it can be extension, file format, such as .txt."
    )


    return parser.parse_args()


def add_replacement_path(file_path, prefix, suffix):
    filepath = Path(file_path)
    filename = filepath.stem  # Get filename without extension
    
    # Extract sample_id from the filename
    sample_id = re.sub(r'_L[0-9]{3}', '', filename)
    sample_id = re.sub(r'_[cs]nv_.*$', '', sample_id)
    sample_id = re.sub(r'_mqc$', '', sample_id)
    

    return prefix + '/' + sample_id + suffix  

    
if __name__ == "__main__":
    # Parse arguments
    args = parse_arguments()

    # Extract arguments
    file_path = args.path
    prefix = args.prefix
    suffix = args.suffix
    
    #print(add_replacement_path(, prefix, suffix))
     