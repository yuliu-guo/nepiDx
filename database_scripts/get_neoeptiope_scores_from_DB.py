'''
Script for calculating neoepitope scores per patient, as number of peptides that pass given thresholds. 
The thresholds are defined in the scores_list file which is a .tsv file with following columns:
  - scoreID - name of the neoepitope score
  - table_name - name of the table based on which the filtering is to be performed (chosen from: v_predictions, v_variant_info_snv, v_variant_info_fusion, v_variant_info_sv)
  - column_name - name of the column to filter on
  - threshold - numeric inclusive threshold for filtering
  - direction - "lower" if peptides should have a score below the threshold, or "higher" if peptides should have a score above the threshold
  - variant_type - which variant type(s) to use (chosen from: "all", "SNV", "fusion", "SV")
Each scoreID can have more than one line, if more than one filter should be applied. Make sure that each line of the same scoreID has the same variant_type.
When filtering using variant tables, use ONLY 1 variant type at a time and make sure the table used matches the variant_type column. 
'''

import sqlite3
import argparse
import pandas as pd
import numpy as np
import os
from functools import reduce

PEPTIDE_TABLE = "v_predictions"

# View mapping
VARIANT_VIEW_MAPPING = {
    "v_variant_info_snv": "SNV",
    "v_variant_info_fusion": "fusion",
    "v_variant_info_sv": "SV"
}

def calculate_neoepitope_scores(con, scores_df):
    """
    Calculates neoepitope scores per subject.
    """
    ### Batch load all required data into memory (Once per unique table)
    tables_cache = {}
        
    for table_name, table_group in scores_df.groupby("table_name"):
        # Identify all unique assessment columns needed for this specific table
        needed_cols = set(table_group["column_name"].unique())
        
        # Always ensure we grab the identifier columns
        if table_name in VARIANT_VIEW_MAPPING:
            needed_cols.update(["subject_id", "code", "treatment_timeline"])
        else:
            needed_cols.update(["subject_id", "code", "treatment_timeline", "variant_type"])
        
        # Wrap column names in double quotes to safely handle special characters like %
        quoted_cols = [f'"{col}"' for col in needed_cols]
        select_clause = ", ".join(quoted_cols)
        
        query = f"""
            SELECT {select_clause}
            FROM {table_name}
            WHERE treatment_timeline = 'pre';
        """
        
        print(f"Loading required columns from table '{table_name}' into memory...")
        df = pd.read_sql_query(query, con)
        
        id_cols = {"subject_id", "code", "treatment_timeline", "variant_type"}
        cols_to_cast = [c for c in df.columns if c not in id_cols]
        df[cols_to_cast] = df[cols_to_cast].apply(pd.to_numeric, errors="coerce")
        
        tables_cache[table_name] = df
        
    if PEPTIDE_TABLE not in tables_cache:
        print(f"Loading the peptides from table '{PEPTIDE_TABLE}' into memory...")
        needed_cols = set(["subject_id", "code", "treatment_timeline", "variant_type"])
        quoted_cols = [f'"{col}"' for col in needed_cols]
        select_clause = ", ".join(quoted_cols)
        
        query = f"""
            SELECT {select_clause}
            FROM {PEPTIDE_TABLE}
            WHERE treatment_timeline = 'pre';
        """
        tables_cache[PEPTIDE_TABLE] = pd.read_sql_query(query , con)

    result_dfs = []
    
    ### Loop through each score
    for score_id, group in scores_df.groupby("scoreID"):
        variant_tables_in_group = [t for t in group["table_name"].unique() if t in VARIANT_VIEW_MAPPING]
        variant_type = group["variant_type"].unique()
        
        if len(variant_type) > 1:
          raise ValueError(f"Make sure the scoreID has the same variant_type value for all rows.") 
          
        if len(variant_tables_in_group) > 1:
          raise ValueError(f"Make sure the score uses up to 1 variant table.") 
        
        if variant_tables_in_group:
            v_table = variant_tables_in_group[0] 
            variant_type_from_table = VARIANT_VIEW_MAPPING[v_table]
            
            if variant_type_from_table != variant_type[0]:
                raise ValueError(f"Variant table does not much the variant type.")
            
            # Fetch cached peptide table dataframe
            pep_df = tables_cache[PEPTIDE_TABLE].copy()

            if pep_df.empty:
                print(f"Warning: Base peptide table missing or empty for scoreID '{score_id}'. Skipping.")
                continue
                
            # Discard all variants of other types immediately from the peptide dataframe
            pep_df = pep_df[pep_df["variant_type"] == variant_type[0]]
            
            # Fetch the target variant view dataframe
            var_df = tables_cache[v_table].copy()
            
            # Perform the relationship join within Pandas on subject_id and code
            working_df = pd.merge(pep_df, var_df, on=["subject_id", "code"], how="inner")
        elif variant_type[0] != "all":
            pep_df = tables_cache[PEPTIDE_TABLE].copy()
            working_df = pep_df[pep_df["variant_type"] == variant_type[0]] 
        else:
            working_df = tables_cache[PEPTIDE_TABLE].copy()
        
        # Start with a boolean mask of all True (same length as the table)
        mask = pd.Series(True, index=working_df.index)
        
        # Chain all condition thresholds together using bitwise AND (&)
        for _, row in group.iterrows():
            col = row["column_name"]
            thresh = row["threshold"]
            direction = str(row["direction"]).strip().lower()
            
            if direction == "lower":
                mask &= (working_df[col] <= thresh)
            elif direction == "higher":
                mask &= (working_df[col] >= thresh)
            else:
                raise ValueError(f"Unknown direction '{direction}' for scoreID '{score_id}'.")
        
        # Apply the combined mask to get rows passing all criteria
        filtered_peptides = working_df[mask]
        
        # Group and count the surviving rows by patient in Python
        if not filtered_peptides.empty:
            df_counts = filtered_peptides.groupby("subject_id").size().reset_index(name=score_id)
        else:
            df_counts = pd.DataFrame(columns=["subject_id", score_id])
            
        result_dfs.append(df_counts)
    
    ### Merge results together 
    # Get all subject_ids
    all_database_subjects = set(tables_cache[PEPTIDE_TABLE]["subject_id"].unique())
        
    final_df = pd.DataFrame({"subject_id": list(all_database_subjects)})
    
    if not result_dfs:
        return final_df
        
    for df_counts in result_dfs:
        final_df = pd.merge(final_df, df_counts, on="subject_id", how="left")
    
    final_df = final_df.fillna(0)
    
    # Convert score columns to integers
    score_cols = final_df.columns.drop("subject_id")
    final_df[score_cols] = final_df[score_cols].astype(int)
    
    return final_df

def main():
  parser = argparse.ArgumentParser(description="Get neoepitope scores per patient")
  parser.add_argument("--database_path", required=True, help="Path to the SQLite database file")
  parser.add_argument("--scores_list", required=True, help="Path to the neoepitope scores list.")
  parser.add_argument("--output_file", required=True, help="Path to the output file with the scores.")

  args = parser.parse_args()
  
  con = sqlite3.connect(args.database_path)

  try:
    # Read in the scores configuration file
    scores_df = pd.read_csv(args.scores_list, sep="\t")
        
    # Call the calculation function
    final_df = calculate_neoepitope_scores(con, scores_df)
        
    if not final_df.empty:
        # Save the final matrix
        final_df.sort_values(by=["subject_id"]).to_csv(args.output_file, sep="\t", index=False)
        print(f"Successfully processed scores and saved results to {args.output_file}")
    else:
        print("Warning: The scores_list file appears to be empty. No scores were processed.")

  except Exception as e:
    print(f"An error occurred while calculating scores: {e}.")

  finally:
    con.close()

if __name__ == "__main__":
    main()