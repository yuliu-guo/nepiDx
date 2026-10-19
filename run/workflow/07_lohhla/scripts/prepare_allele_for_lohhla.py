#!/usr/bin/env python3

import argparse
import re
from pathlib import Path
from typing import Optional


def load_p_group_map(p_group_file: Path) -> dict:
    """
    Parse an HLA P-group file.

    Expected format, semicolon-separated:

        A*;01:01:01:01/01:01:01:02;01:01P

    Returns mapping like:

        A*01:01:01:01 -> A*01:01
        A*01:01:01:02 -> A*01:01
    """
    allele_to_pgroup = {}

    with p_group_file.open("r") as infile:
        for line in infile:
            line = line.strip()

            if not line or line.startswith("#"):
                continue

            fields = line.split(";")

            if len(fields) < 3:
                continue

            gene = fields[0].strip()          # e.g. A*
            alleles = fields[1].split("/")    # e.g. 02:01:01:01/02:01:01:02
            p_group = fields[2].strip()       # e.g. 02:01P

            if not gene or not p_group:
                continue

            if p_group.endswith("P"):
                p_group_root = p_group[:-1]
            else:
                p_group_root = p_group

            for allele in alleles:
                allele = allele.strip()

                if not allele:
                    continue

                allele_to_pgroup[f"{gene}{allele}"] = f"{gene}{p_group_root}"

    return allele_to_pgroup


def map_to_p_group(full_allele: str, pgroup_map: dict) -> str:
    """
    Map:

        HLA-A*02:01:01:01

    to its P-group allele if available, for example:

        HLA-A*02:01

    If no P-group mapping is found, return the original allele.
    """
    match = re.match(r"^(HLA-([ABC])\*)(\S+)$", full_allele)

    if not match:
        return full_allele

    full_prefix = match.group(1)       # HLA-A*
    locus = match.group(2)             # A
    allele_part = match.group(3)       # 02:01:01:01

    pgroup_key = f"{locus}*{allele_part}"  # A*02:01:01:01

    if pgroup_key in pgroup_map:
        return f"HLA-{pgroup_map[pgroup_key]}"

    return full_allele


def trim_to_two_field_resolution(full_allele: str) -> Optional[str]:
    """
    Convert:

        HLA-A*02:01:01:01
        HLA-A*02:01

    to:

        hla_a_02_01

    Only HLA-A, HLA-B, and HLA-C are retained.
    """
    match = re.match(r"^HLA-([ABC])\*(\S+)$", full_allele)

    if not match:
        return None

    locus = match.group(1).lower()
    allele_part = match.group(2)

    allele_fields = allele_part.split(":")

    if len(allele_fields) < 2:
        return None

    first_two_fields = allele_fields[:2]

    return f"hla_{locus}_{first_two_fields[0]}_{first_two_fields[1]}"


def convert_allele(allele: str, pgroup_map: dict) -> Optional[str]:
    """
    Full conversion:

        HLA-A*02:01:01:01
        -> map to P-group
        -> trim to two-field resolution
        -> hla_a_02_01
    """
    mapped_allele = map_to_p_group(allele, pgroup_map)
    converted = trim_to_two_field_resolution(mapped_allele)

    return converted


def main():
    parser = argparse.ArgumentParser(
        description="Prepare HLA allele TSV for LOHHLA using P-group mapping."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input allele TSV file."
    )
    parser.add_argument(
        "--p-group-file",
        required=True,
        help="HLA P-group mapping file."
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output LOHHLA-formatted allele TSV file."
    )
    parser.add_argument(
        "--quality-threshold",
        type=int,
        default=0,
        help="Minimum HLA allele quality score to keep. Default: 0."
    )

    args = parser.parse_args()

    input_path = Path(args.input)
    p_group_path = Path(args.p_group_file)
    output_path = Path(args.output)

    output_path.parent.mkdir(parents=True, exist_ok=True)

    pgroup_map = load_p_group_map(p_group_path)

    converted_alleles = []
    seen = set()

    hla_pattern = re.compile(r"^(HLA-[ABC]\*\S+)\s+(\d+)$")

    with input_path.open("r") as infile:
        for line in infile:
            line = line.strip()

            if not line:
                continue

            match = hla_pattern.match(line)

            if not match:
                continue

            allele = match.group(1)
            quality = int(match.group(2))

            if quality < args.quality_threshold:
                continue

            converted = convert_allele(allele, pgroup_map)

            if converted is None:
                continue

            if converted not in seen:
                converted_alleles.append(converted)
                seen.add(converted)

    with output_path.open("w") as outfile:
        for allele in converted_alleles:
            outfile.write(f"{allele}\n")


if __name__ == "__main__":
    main()