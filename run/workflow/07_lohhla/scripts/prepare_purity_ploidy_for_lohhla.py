#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(
        description="Prepare PURPLE purity/ploidy file for LOHHLA."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input PURPLE purity TSV file."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output LOHHLA-formatted purity/ploidy TSV file."
    )
    parser.add_argument(
        "--sample-name",
        required=True,
        help="Sample name to use as the first column in the LOHHLA file."
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with input_path.open("r", newline="") as infile:
        reader = csv.DictReader(infile, delimiter="\t")
        rows = list(reader)

    if len(rows) != 1:
        raise ValueError(
            f"Expected exactly one data row in {input_path}, found {len(rows)}"
        )

    row = rows[0]

    tumor_purity = row["purity"]
    tumor_ploidy = row["ploidy"]

    with output_path.open("w") as outfile:
        outfile.write("\tPloidy\ttumorPurity\ttumorPloidy\n")
        outfile.write(f"{args.sample_name}\t2\t{tumor_purity}\t{tumor_ploidy}\n")


if __name__ == "__main__":
    main()
