import os
import subprocess
import argparse
import csv
import re
import csv


def create_mhchammer_inventory(patient_id, normal_wxs, tumor_wxs, tumor_rna, purple_res, scratch_dir):
    """
    Creates MHC Hammer inventory file for a given patient.

    Parameters:
        patient_id (str): The patient identifier
        normal_wxs (str): Path to the normal WXS BAM file
        tumor_wxs (str): Path to the tumor WXS BAM file
        tumor_rna (str): Path to the tumor RNA BAM file
        scratch_dir (str): Directory to save the inventory file
    """
    ## Get ploidy and purity
    with open(purple_res, mode='r', newline='') as f:
      reader = csv.DictReader(f, delimiter='\t')
      first_row = next(reader)
      
      ploidy = first_row['ploidy']
      purity = first_row['purity']
    
    inventory_path = os.path.join(scratch_dir, f"{patient_id}_mhchammer_inventory.csv")
    #allele_path = os.path.join(scratch_dir, f"{patient_id}_alleles.csv")
    with open(inventory_path, 'w+') as f:
        writer = csv.writer(f, delimiter=',')
        writer.writerow(["patient", "sample_name", "sample_type", "bam_path", "sequencing_type", "purity", "ploidy", "normal_sample_name"])
        writer.writerow([patient_id, f"{patient_id}_normal_wxs", "normal", normal_wxs, "wxs", purity, ploidy, f"{patient_id}_normal_wxs"])
        writer.writerow([patient_id, f"{patient_id}_tumor_wxs", "tumor", tumor_wxs, "wxs", purity, ploidy, f"{patient_id}_normal_wxs"])
        writer.writerow([patient_id, f"{patient_id}_tumor_rna", "tumor", tumor_rna, "rna", purity, ploidy, f"{patient_id}_normal_wxs"])
    return inventory_path

def main():
    parser = argparse.ArgumentParser(description="Run MHCflurry predictions")
    parser.add_argument("-p", "--patient_id", required=True, help="Patient ID")
    parser.add_argument("--normal_wxs", required=True, help="Path to normal WXS BAM file")
    parser.add_argument("--tumor_wxs", required=True, help="Path to tumor WXS BAM file")
    parser.add_argument("--tumor_rna", required=True, help="Path to tumor RNA BAM file")
    parser.add_argument("--purple_result_file", required=True, help="Path to file with purity and ploidy estimates from PURPLE.")
    #parser.add_argument("--hla_alleles", required=True, help="Path to HLA allele predictions")
    parser.add_argument("--scratch_dir", type=str, required=True, help="Path to the scratch directory")
    parser.add_argument("-q", "--quality_threshold", type=int, default=0, help="Quality threshold for filtering alleles")
    
    args = parser.parse_args()
    
    #extract_and_transform_alleles(
    #  args.patient_id,
    #  args.hla_alleles, 
    #  args.scratch_dir, 
    #  args.quality_threshold
    #)
    
    inventory_file = create_mhchammer_inventory(
        args.patient_id,
        args.normal_wxs,
        args.tumor_wxs,
        args.tumor_rna,
        args.purple_result_file,
        args.scratch_dir
    )
    print(f"Inventory file created at: {inventory_file}")

if __name__ == "__main__":
    main()