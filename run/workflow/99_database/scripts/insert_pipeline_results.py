import sqlite3
import pandas as pd
import argparse
import datetime
import os
import re
import csv

def load_p_group_map(p_group_file):
    """
    Parses the P group file and returns a mapping of full allele (e.g. 01:315)
    to its P group root (e.g. 01:03), using gene prefix (e.g. A*).
    Skips entries with missing P-group root values.
    """
    allele_to_pgroup = {}

    with open(p_group_file, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue

            fields = line.strip().split(';')
            if len(fields) < 3:
                continue

            gene = fields[0].strip()           # e.g., 'A*'
            alleles = fields[1].split('/')     # e.g., ['01:315', '01:01:153']
            p_group = fields[2].strip()        # e.g., '01:03P'

            if not p_group:                    # Skip if no P group root
                continue

            # Extract root from p_group (removing trailing 'P')
            root = p_group.replace('P', '')

            for a in alleles:
                trimmed = a.strip()
                if trimmed:
                    allele_to_pgroup[f'{gene}{trimmed}'] = f'{gene}{root}'

    return allele_to_pgroup


def map_to_p_group(full_allele, pgroup_map):
    """
    Maps full allele to its HLA P-group root (e.g., HLA-A*01:315 -> HLA-A*01:03),
    or returns the original allele if no mapping exists.
    """
    match = re.match(r'^(HLA-[A-Z]\*)(\S+)', full_allele)
    if not match:
      return full_allele  # fallback - return full allele even if pattern is weird

    prefix, allele_part = match.groups()       # e.g., 'HLA-A*', '01:12'
    pgroup_key = f'{prefix[4:]}{allele_part}'  # strip 'HLA-' -> 'A*01:12'
    
    if pgroup_key in pgroup_map:
      return f'HLA-{pgroup_map[pgroup_key]}'
    else:
      return full_allele                      # no mapping found, return as-is

def get_variant_type_from_filename(filename):
    """
    Extracts the variant type from a given filename.

    Parameters:
        filename (str): The filename to extract the variant type from
    Returns:
        str: The extracted variant type
    """
    if "snv" in filename.lower():
        return "SNV"
    elif "fusion" in filename.lower():
        return "fusion"
    else:
        return "SV"

def find_dataset_id(con, subject_id, glass_nr):
    """
    Finds the dataset ID for a given subject ID and glass number.

    Parameters:
        con (sqlite3.Connection): SQLite database connection
        subject_id (str): The subject identifier
        glass_nr (str): The glass number
    Returns:
        int: The dataset ID if found, None otherwise
    """
    cur = con.cursor()
    dataset_search = cur.execute('SELECT id FROM dataset WHERE subject_id = ? AND glass_nr = ?', (subject_id, glass_nr)).fetchone()
    if dataset_search is None:
        raise ValueError(f"No dataset found for subject_id '{subject_id}' and glass_nr '{glass_nr}'. Please ensure the dataset exists in the database before inserting results.")   
    else:
        return dataset_search[0]

### Insert into allele table ###
def insert_allele_typing_results(con, results_dir, current_datasets, dataset_id, p_group_file, with_update=False):
    """
    Inserts HLA allele typing results into the database for a given patient.

    Parameters:
        con (sqlite3.Connection): SQLite database connection
        results_dir (str): Directory containing the pipeline results files
        current_datasets (set): Set of dataset IDs already present in the database
        dataset_id (int): The dataset ID for which to insert allele typing results
        with_update (bool): Whether to update existing allele typing results if the dataset already exist in the database (default: False)
    """
    # Get filepath for allele typing results
    allele_typing_dir = os.path.join(results_dir, "01_hla-typing")
    for file in os.listdir(allele_typing_dir):
        if file.endswith("_tumor_DNA_allele.tsv"):
            allele_typing_file = os.path.join(allele_typing_dir, file)
            break
    
    allele_typing_df = pd.read_csv(allele_typing_file, sep=' ', header=None, names=['allele', 'quality'])
    p_group_map = load_p_group_map(p_group_file)
    cur = con.cursor()

    # Insert allele typing results into the database
    print(f"Starting insert of allele typing results for dataset_id: {dataset_id}")
    if with_update and dataset_id in current_datasets:
        print(f"Allele typing results for dataset_id '{dataset_id}' already exist in the database. Removing existing records.")
        cur.execute('''
            DELETE FROM allele
            WHERE dataset_id = ?
        ''', (dataset_id,))
    for _, row in allele_typing_df.iterrows():
        if dataset_id not in current_datasets or with_update:
            allele_full_name = row['allele']
            allele_name = map_to_p_group(allele_full_name, p_group_map)
            cur.execute('''
            INSERT INTO allele (dataset_id, allele_full_name, allele_name, quality, created_at)
            VALUES (?, ?, ?, ?, ?)
            ''', (dataset_id, allele_full_name, allele_name, row['quality'], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        else:
            print(f"Allele typing results for dataset_id '{dataset_id}' already exist in the database. Skipping insertion.")
            break
    print(f"Finished insert of allele typing results for dataset_id: {dataset_id}")
    
### Insert into allele_loh table ###
def insert_allele_loh_results(con, results_dir, current_datasets, dataset_id, with_update=False):
    """
    Inserts HLA allele LOH results into the database for a given patient.

    Parameters:
        con (sqlite3.Connection): SQLite database connection
        results_dir (str): Directory containing the pipeline results files
        current_datasets (set): Set of dataset IDs already present in the database
        dataset_id (int): The dataset ID for which to insert allele typing results
        with_update (bool): Whether to update existing allele typing results if the dataset already exist in the database (default: False)
    """
    # Get filepath for allele typing results
    allele_loh_file = os.path.join(results_dir, "07_spechla/merge.hla.copy.txt")
    
    if not os.path.isfile(allele_loh_file):
        print("HLA LOH results not available for the patient. Skipping insert.")
        return
    
    allele_loh_df = pd.read_csv(allele_loh_file, sep='\t')
    cur = con.cursor()

    # Insert allele loh results into the database
    print(f"Starting insert of allele LOH results for dataset_id: {dataset_id}")
    if with_update and dataset_id in current_datasets:
        print(f"Allele LOH results for dataset_id '{dataset_id}' already exist in the database. Removing existing records.")
        cur.execute('''
            DELETE FROM allele_loh
            WHERE dataset_id = ?
        ''', (dataset_id,))
        
    for _, row in allele_loh_df.iterrows():
        if dataset_id not in current_datasets or with_update:
            is_loh = "yes" if row["LOH"] == "Y" else "no"
            
            cur.execute('''
            INSERT INTO allele_loh (dataset_id, HLA, allele_name_1, allele_name_2, copyratio, 
                                    kept_HLA, loss_HLA, allele_freq_1, allele_freq_2, het_num, is_LOH, created_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset_id, row['HLA'], row['Allele1'], row['Allele2'], row['copyratio'], row['KeptHLA'], 
                  row['LossHLA'], row['Freq1'], row['Freq2'], row['Het_num'], is_loh, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        else:
            print(f"Allele LOH results for dataset_id '{dataset_id}' already exist in the database. Skipping insertion.")
            break
    print(f"Finished insert of allele LOH results for dataset_id: {dataset_id}")

def insert_variant_info(con, results_dir, current_datasets, dataset_id, with_update=False):
    """
    Inserts variant information into the database for a given patient.

    Parameters:
        con (sqlite3.Connection): SQLite database connection
        results_dir (str): Directory containing the pipeline results files
        current_datasets (set): Set of dataset IDs already present in the database
        dataset_id (int): The dataset ID for which to insert variant information
        with_update (bool): Whether to update existing variant information if the dataset already exist in the database (default: False)
    """
    gene_col_dict = {
        "SNV": "SYMBOL",
        "fusion": "gene_name",
        "SV": "sv_gene"
    }
    sql_col_dict = {
        "SNV": "snv_variant_number",
        "fusion": "fusion_variant_number",
        "SV": "sv_variant_number"
    }
    cur = con.cursor()
    
    print(f"Inserting variant summary results and variant number for dataset_id: {dataset_id}")
    # Get filepaths for variant summary results
    files = os.listdir(os.path.join(results_dir, "00_variants_summary"))
    for file in files:
        if dataset_id not in current_datasets or with_update:
            path = os.path.join(results_dir, "00_variants_summary", file)
            variant_type = get_variant_type_from_filename(file)
            variant_summary_df = pd.read_csv(path, sep='\t')
            
            # Insert variant number into dataset table
            variant_number = len(variant_summary_df)
            cur.execute(f'''
                UPDATE dataset
                SET {sql_col_dict[variant_type]} = ?
                WHERE id = ?
            ''', (variant_number, dataset_id))
                
            for _, row in variant_summary_df.iterrows():
                gene = row[gene_col_dict[variant_type]]
                code = row['ID']
                # Insert into variant table
                cur.execute('''
                    INSERT INTO variant (code, variant_type, dataset_id, gene, created_at)
                    VALUES (?, ?, ?, ?, ?)
                ''', (code, variant_type, dataset_id, gene, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
                
                # Insert into variant_info table
                variant_id = cur.lastrowid
                for col_name, col_value in row.items():
                    if col_name not in ['ID', gene_col_dict[variant_type]]:
                        cur.execute('''
                            INSERT INTO variant_info (variant_id, attribute_name, attribute_value, created_at)
                            VALUES (?, ?, ?, ?)
                        ''', (variant_id, col_name, str(col_value), datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        else:
            print(f"Variants for dataset_id '{dataset_id}' already exist in the database. Skipping insertion.")
            break
                    
    print(f"Finished insert of variant summary results for dataset_id: {dataset_id}")
    
def insert_mutated_peptide_results(con, results_dir, current_datasets, dataset_id, with_update=False):
    """
    Inserts peptides into the database for a given patient.

    Parameters:
        con (sqlite3.Connection): SQLite database connection
        results_dir (str): Directory containing the pipeline results files
        current_datasets (set): Set of dataset IDs already present in the database
        dataset_id (int): The dataset ID for which to insert peptide prediction results
        with_update (bool): Whether to update existing peptide prediction results if the dataset already exist in the database (default: False)
    """
    cur = con.cursor()
    
    sql_col_dict = {
        "SNV": "snv_peptide_number",
        "fusion": "fusion_peptide_number",
        "SV": "sv_peptide_number"
    }
    
    # Mapping variants: code -> id
    variants = {row[0]: row[1] for row in cur.execute('SELECT code, id FROM variant WHERE dataset_id = ?', (dataset_id,)).fetchall()}
    
    print(f"Inserting peptides and peptide number for dataset_id: {dataset_id}")
    peptide_dir = os.path.join(results_dir, "02_mutated_peptide")
    files = os.listdir(peptide_dir)
    
    for file in files:
        if file.endswith("_similarity.csv"):
            path = os.path.join(peptide_dir, file)
            peptide_df = pd.read_csv(path, sep='\t')
            peptide_df["ID"] = peptide_df["ID"].str.split("_").str[0].str.split("|").str[0]
            peptide_df = peptide_df.drop_duplicates(subset=["ID", "Mut_peptide"], keep="first")
            variant_type = get_variant_type_from_filename(file)
            
            # Insert peptide number into dataset table
            peptide_number = len(peptide_df)
            cur.execute(f'''
                UPDATE dataset
                SET {sql_col_dict[variant_type]} = ?
                WHERE id = ?
            ''', (peptide_number, dataset_id))
            
            for _, row in peptide_df.iterrows():
                variant_code = row['ID']                  
                variant_id = variants.get(variant_code)
                if not variant_id:
                    raise ValueError(f"No variant found in the database with code '{variant_code}' for dataset_id '{dataset_id}'. Please ensure the variant information is inserted into the database before inserting peptide prediction results.")
                
                if dataset_id not in current_datasets or with_update:
                    # Insert into peptide table
                    cur.execute('''
                        INSERT INTO peptide (sequence, wildtype_sequence, variant_id, hamming_similarity, blosum62_2_similarity, created_at)
                        VALUES (?, ?, ?, ?, ?, ?)
                    ''', (row['Mut_peptide'], row['WT_peptide'], variant_id, row['hamming_similarity'], row['blosum62_2_similarity'], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
                else:
                    print(f"Peptides for dataset_id '{dataset_id}' already exist in the database. Skipping insertion of file {file}.")
                    break
            
    print(f"Finished insert of peptide prediction results for dataset_id: {dataset_id}")

def insert_peptide_prediction_results(con, results_dir, current_datasets, dataset_id, with_update=False):
    """
    Inserts peptide binding and immunogenicity prediction results into the database for a given patient.

    Parameters:
        con (sqlite3.Connection): SQLite database connection
        results_dir (str): Directory containing the pipeline results files
        current_datasets (set): Set of dataset IDs already present in the database
        dataset_id (int): The dataset ID for which to insert peptide prediction results
        with_update (bool): Whether to update existing peptide prediction results if the dataset already exist in the database (default: False)
    """
    cur = con.cursor()
    
    # Mapping variants: code -> id
    variants = {row[0]: row[1] for row in cur.execute('SELECT code, id FROM variant WHERE dataset_id = ?', (dataset_id,)).fetchall()}
    # Mapping peptides: (variant_id, sequence) -> id
    peptide_query = f'''
        SELECT variant_id, sequence, id FROM peptide 
        WHERE variant_id IN ({",".join(["?"] * len(variants))})
    '''
    peptides = {(row[0], row[1]): row[2] for row in cur.execute(peptide_query, list(variants.values())).fetchall()}
    # Mapping alleles: name -> id
    alleles = {row[0]: row[1] for row in cur.execute('SELECT allele_name, id FROM allele WHERE dataset_id = ?', (dataset_id,)).fetchall()}
    
    print(f"Inserting peptide binding and immunogenicity prediction results for dataset_id: {dataset_id}")
    peptide_pred_dir = os.path.join(results_dir, "03_pMHC_prediction")
    immunogenicity_dir = os.path.join(results_dir, "08_immunogenicity")
    
    prediction_files = os.listdir(peptide_pred_dir)
    immunogenicity_files = os.listdir(immunogenicity_dir)
    
    paths_to_insert = []
    
    for file in prediction_files:
        is_from_variant = file.startswith("fusion") or file.startswith("snv") or file.startswith("sv")
        if is_from_variant and file.endswith("_with_DAI.csv"):
            paths_to_insert.append(os.path.join(peptide_pred_dir, file))
            
    for file in immunogenicity_files:
        is_from_variant = file.startswith("fusion") or file.startswith("snv") or file.startswith("sv")
        if is_from_variant:
            paths_to_insert.append(os.path.join(immunogenicity_dir, file))
    
    for path in paths_to_insert:
        print("Inserting file: ", path)
        prediction_df = pd.read_csv(path, sep='\t')
        prediction_df["Identity"] = prediction_df["Identity"].str.split("_").str[0].str.split("|").str[0]
        prediction_df = prediction_df.drop_duplicates(subset=["HLA_allele", "Peptide", "Identity"], keep="first")
        
        for _, row in prediction_df.iterrows():
            variant_code = row['Identity']            
            variant_id = variants.get(variant_code)
            if not variant_id:
                raise ValueError(f"No variant found in the database with code '{variant_code}' for dataset_id '{dataset_id}'. Please ensure the variant information is inserted into the database before inserting peptide prediction results.")
            
            if dataset_id not in current_datasets or with_update:
                # Insert into peptide2prediction table 
                peptide_id = peptides.get((variant_id, row['Peptide']))
                if not peptide_id:
                    raise ValueError(f"No peptide found in the database with variant_id '{variant_id}' and sequence '{row['Peptide']}'. Please ensure the peptide information is inserted into the database before inserting peptide prediction results.")
                
                allele_name = row['HLA_allele']
                allele_id = alleles.get(allele_name)
                if not allele_id:
                    allele_name = allele_name[0:5] + '*' + allele_name[5:]
                    allele_id = alleles.get(allele_name)
                    if not allele_id:
                        raise ValueError(f"No allele found in the database with name '{allele_name}' for dataset_id '{dataset_id}'. Please ensure the allele typing results are inserted into the database before inserting peptide prediction results.")
                    
                prediction_tool = row.get('Prediction_tool', row.get('Immunogenicity_tool', 'netMHCstabpan')) # There was a bug in netMHCstabpan DAI code, now it's fixed; change the 'netMHCstabpan' to 'Unknown' after the pipeline has been re-run
                for col_name, col_value in row.items():
                    if col_name not in ['Patient_ID', 'Identity', 'Peptide', 'HLA_allele', 'Prediction_tool', 'Immunogenicity_tool']:
                        cur.execute('''
                            INSERT INTO peptide2prediction (peptide_id, allele_id, prediction_tool, prediction_type, value, created_at)
                            VALUES (?, ?, ?, ?, ?, ?)
                        ''', (peptide_id, allele_id, prediction_tool, col_name, col_value, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            else:
                print(f"Peptide prediction results for dataset_id '{dataset_id}' already exist in the database. Skipping insertion of file {file}.")
                break
            
    print(f"Finished insert of peptide binding and immunogenicity prediction results for dataset_id: {dataset_id}")

def insert_variant_and_peptide_results(con, results_dir, current_datasets, dataset_id, with_update=False):
    """
    Inserts variant summary results into the database for a given patient.

    Parameters:
        con (sqlite3.Connection): SQLite database connection
        results_dir (str): Directory containing the pipeline results files
        current_datasets (set): Set of dataset IDs already present in the database
        dataset_id (int): The dataset ID for which to insert variant summary results
        with_update (bool): Whether to update existing variant summary results if the dataset already exist in the database (default: False)
    """   
    cur = con.cursor()
    
    print(f"Starting insert of variant and peptide results for dataset_id: {dataset_id}")
    if with_update and dataset_id in current_datasets:
        print(f"Variant and peptide results for dataset_id '{dataset_id}' already exist in the database. Removing existing records.")
        
        print(f"Removing existing records from 'peptide2prediction' table for dataset_id '{dataset_id}'")
        cur.execute('''
            DELETE FROM peptide2prediction
            WHERE peptide_id IN (SELECT id FROM peptide WHERE variant_id IN (SELECT id FROM variant WHERE dataset_id = ?))
        ''', (dataset_id,))
        
        print(f"Removing existing records from 'peptide' table for dataset_id '{dataset_id}'")
        cur.execute('''
            DELETE FROM peptide
            WHERE variant_id IN (SELECT id FROM variant WHERE dataset_id = ?)
        ''', (dataset_id,))
        
        print(f"Removing existing records from 'variant_info' table for dataset_id '{dataset_id}'")
        cur.execute('''
            DELETE FROM variant_info
            WHERE variant_id IN (SELECT id FROM variant WHERE dataset_id = ?)
        ''', (dataset_id,))
        
        print(f"Removing existing records from 'variant' table for dataset_id '{dataset_id}'")
        cur.execute('''
            DELETE FROM variant
            WHERE dataset_id = ?
        ''', (dataset_id,))
        
    ### Insert into variant and variant_info tables ###
    insert_variant_info(con, results_dir, current_datasets, dataset_id, with_update=with_update)
    
    ### Insert into peptide table ###
    insert_mutated_peptide_results(con, results_dir, current_datasets, dataset_id, with_update=with_update)
    
    ### Insert into peptide2prediction table ###
    insert_peptide_prediction_results(con, results_dir, current_datasets, dataset_id, with_update=with_update)
    
    print(f"Finished insert of variant and peptide results for dataset_id: {dataset_id}")
    
### Insert into peptide_ms table
def insert_peptide_ms_results(con, results_dir, current_datasets, dataset_id, with_update=False):
    """
    Inserts DIA-NN results into the database for a given patient.

    Parameters:
        con (sqlite3.Connection): SQLite database connection
        results_dir (str): Directory containing the pipeline results files
        current_datasets (set): Set of dataset IDs already present in the database
        dataset_id (int): The dataset ID for which to insert variant summary results
        with_update (bool): Whether to update existing variant summary results if the dataset already exist in the database (default: False)
    """ 
    # Get filepath for summarized MS peptide results
    peptide_ms_dir = os.path.join(results_dir, "09_diann")
    
    if not os.path.isdir(peptide_ms_dir):
        print("Results for peptides from MS are not available for the patient. Skipping insert.")
        return
    
    ms_file_list = os.listdir(peptide_ms_dir)
    ms_file = next((f for f in ms_file_list if f.endswith("_summary_specific_peptides.tsv")), None)
    
    if ms_file is None:
        print("Summary file is missing from the DIA-NN results. Make sure it was generated correctly. Skipping insert.")
        return
    
    ms_filepath = os.path.join(peptide_ms_dir, ms_file)
    summary_ms_df = pd.read_csv(ms_filepath, sep='\t')
    cur = con.cursor()
    
    # Handle update - delete previous results
    if with_update and dataset_id in current_datasets:
        print(f"Results for peptides from MS for dataset_id '{dataset_id}' already exist in the database. Removing existing records.")
        cur.execute('''
            DELETE FROM peptide_ms
            WHERE dataset_id = ?
        ''', (dataset_id,))
        
    # Insert results into the database
    print(f"Starting insert of results for peptides from MS for dataset_id: {dataset_id}")
    for _, row in summary_ms_df.iterrows():
        if dataset_id not in current_datasets or with_update:
            cur.execute('''
            INSERT INTO peptide_ms (dataset_id, code, sequence, modified_sequence, 
                                    n_runs, n_precursors, protein_ids, protein_names, 
                                    genes, precursor_norm_max, precursor_norm_mean,
                                    precursor_norm_median, n_proteins_matched, peptide_type,
                                    created_at)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (dataset_id, row['id'], row['sequence'], row['modified_sequence'], 
                  row['n_runs'], row['n_precursors'], row['protein_ids'], row['protein_names'], 
                  row['genes'], row['precursor_norm_max'], row['precursor_norm_mean'], 
                  row['precursor_norm_median'], row['n_proteins_matched'], row['type'], 
                  datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        else:
            print(f"Results for peptides from MS for dataset_id '{dataset_id}' already exist in the database. Skipping insertion.")
            break
    print(f"Finished insert of the results for peptides from MS for dataset_id: {dataset_id}")
    
    
def insert_peptide_ms_prediction_results(con, results_dir, current_datasets, dataset_id, with_update=False):
    """
    Inserts peptide from MS binding and immunogenicity prediction results into the database for a given patient.

    Parameters:
        con (sqlite3.Connection): SQLite database connection
        results_dir (str): Directory containing the pipeline results files
        current_datasets (set): Set of dataset IDs already present in the database
        dataset_id (int): The dataset ID for which to insert peptide prediction results
        with_update (bool): Whether to update existing peptide prediction results if the dataset already exist in the database (default: False)
    """
    cur = con.cursor()
    
    # Mapping MS peptides: code -> id
    peptides = {row[0]: row[1] for row in cur.execute('SELECT code, id FROM peptide_ms WHERE dataset_id = ?', (dataset_id,)).fetchall()}

    # Mapping alleles: name -> id
    alleles = {row[0]: row[1] for row in cur.execute('SELECT allele_name, id FROM allele WHERE dataset_id = ?', (dataset_id,)).fetchall()}
    
    peptide_pred_dir = os.path.join(results_dir, "03_pMHC_prediction")
    immunogenicity_dir = os.path.join(results_dir, "08_immunogenicity")
    
    prediction_files = os.listdir(peptide_pred_dir)
    immunogenicity_files = os.listdir(immunogenicity_dir)
    
    paths_to_insert = []
    
    for file in prediction_files:
        is_from_ms = file.startswith("ms8mer") or file.startswith("ms9mer")
        if is_from_ms:
            paths_to_insert.append(os.path.join(peptide_pred_dir, file))
            
    for file in immunogenicity_files:
        is_from_ms = file.startswith("ms8mer") or file.startswith("ms9mer")
        if is_from_ms:
            paths_to_insert.append(os.path.join(immunogenicity_dir, file))
    
    if len(paths_to_insert) == 0:
        print("No MS peptide prediction results available for dataset_id: {dataset_id}. Skipping insert.")
        return
        
    print(f"Inserting MS peptide binding and immunogenicity prediction results for dataset_id: {dataset_id}")
    
    if with_update and dataset_id in current_datasets:
        print(f"MS peptide prediction results for dataset_id '{dataset_id}' already exist in the database. Removing existing records.")
        
        cur.execute('''
            DELETE FROM peptide_ms2prediction
            WHERE peptide_id IN (SELECT id FROM peptide_ms WHERE dataset_id = ?)
        ''', (dataset_id,))
    
    for path in paths_to_insert:
        print("Inserting file: ", path)
        prediction_df = pd.read_csv(path, sep='\t')
        prediction_df["Identity"] = prediction_df["Identity"].str.split("_").str[0].str.split("|").str[0]
        prediction_df = prediction_df.drop_duplicates(subset=["HLA_allele", "Peptide", "Identity"], keep="first")
        
        for _, row in prediction_df.iterrows():
            peptide_code = row['Identity']            
            peptide_id = peptides.get(peptide_code)
            if not peptide_id:
                raise ValueError(f"No MS peptide found in the database with code '{peptide_code}' for dataset_id '{dataset_id}'. Please ensure the MS peptide information is inserted into the database before inserting MS peptide prediction results.")
            
            if dataset_id not in current_datasets or with_update:
                # Insert into peptide_ms2prediction table 
                allele_name = row['HLA_allele']
                allele_id = alleles.get(allele_name)
                if not allele_id:
                    allele_name = allele_name[0:5] + '*' + allele_name[5:]
                    allele_id = alleles.get(allele_name)
                    if not allele_id:
                        raise ValueError(f"No allele found in the database with name '{allele_name}' for dataset_id '{dataset_id}'. Please ensure the allele typing results are inserted into the database before inserting peptide prediction results.")
                    
                prediction_tool = row.get('Prediction_tool', row.get('Immunogenicity_tool', 'netMHCstabpan')) # There was a bug in netMHCstabpan DAI code, now it's fixed; change the 'netMHCstabpan' to 'Unknown' after the pipeline has been re-run
                
                nmer_sequence = row['Peptide']
                
                for col_name, col_value in row.items():
                    if col_name not in ['Patient_ID', 'Identity', 'Peptide', 'HLA_allele', 'Prediction_tool', 'Immunogenicity_tool']:
                        cur.execute('''
                            INSERT INTO peptide_ms2prediction (peptide_id, nmer_sequence, allele_id, prediction_tool, prediction_type, value, created_at)
                            VALUES (?, ?, ?, ?, ?, ?, ?)
                        ''', (peptide_id, nmer_sequence, allele_id, prediction_tool, col_name, col_value, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            else:
                print(f"MS peptide prediction results for dataset_id '{dataset_id}' already exist in the database. Skipping insertion of file {file}.")
                break
            
    print(f"Finished insert of MS peptide binding and immunogenicity prediction results for dataset_id: {dataset_id}")
                        
### Insert into gene_signature_value table ###
def insert_gene_signature_values(con, results_dir, current_datasets, dataset_id, with_update=False):
    """
    Inserts gene signature values into the database for a given patient.
    Parameters:
        con (sqlite3.Connection): SQLite database connection
        results_dir (str): Directory containing the pipeline results files
        current_datasets (set): Set of dataset IDs already present in the database
        dataset_id (int): The dataset ID for which to insert gene signature values
        with_update (bool): Whether to update existing gene signature values if the dataset already exist in the database (default: False)
    """
    gene_sig_dir = os.path.join(results_dir, "05_TME")
    for filepath in os.listdir(gene_sig_dir):
      if "_gene_signatures.csv" in filepath:
        gene_sig_filepath = filepath
        break

    gene_signatures_df = pd.read_csv(os.path.join(gene_sig_dir, gene_sig_filepath), sep='\t', names=['score', 'method', 'value'])
    gene_signatures_df = gene_signatures_df.drop([0])
    gene_signatures_df = gene_signatures_df[['score', 'value']]
    cur = con.cursor()
    
    print(f"Starting insert of gene signature values for dataset_id: {dataset_id}")
    if with_update and dataset_id in current_datasets:
        print(f"Gene signature values for dataset_id '{dataset_id}' already exist in the database. Removing existing records.")
        cur.execute('''
            DELETE FROM gene_signature_value
            WHERE dataset_id = ?
        ''', (dataset_id,))
    
    for _, row in gene_signatures_df.iterrows():
        if dataset_id not in current_datasets or with_update:
            gene_signature_type = cur.execute('SELECT name FROM gene_signature_type WHERE name = ?', (row['score'],)).fetchone()
            if gene_signature_type is None:
                raise ValueError(f"No gene signature type found in the database with name '{row['score']}'. Please ensure the gene signature types are inserted into the database before inserting gene signature values.")
            cur.execute('''
                INSERT INTO gene_signature_value (dataset_id, gene_signature_type, value, created_at)
                VALUES (?, ?, ?, ?)
            ''', (dataset_id, row['score'], row['value'], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        else:
            print(f"Gene signature values for dataset_id '{dataset_id}' already exist in the database. Skipping insertion.")
            break
    print(f"Finished insert of gene signature values for dataset_id: {dataset_id}")

### Insert ploidy and purity into dataset table ###  
def insert_ploidy_and_purity(con, purple_results, current_datasets, dataset_id, with_update=False):
    # Get ploidy and purity
    with open(purple_results, mode='r', newline='') as f:
      reader = csv.DictReader(f, delimiter='\t')
      first_row = next(reader)
      
      ploidy = first_row['ploidy']
      purity = first_row['purity']
      
    cur = con.cursor()
        
    if dataset_id not in current_datasets or with_update:
        if with_update and dataset_id in current_datasets:
            print(f"Ploidy and purity for dataset_id '{dataset_id}' already exist in the database. Updating existing records.")
            
        cur.execute('''
            UPDATE dataset
            SET purity = ?, ploidy = ?
            WHERE id = ?
        ''', (purity, ploidy, dataset_id))
    else:
        print(f"Ploidy and purity for dataset_id '{dataset_id}' already exist in the database. Skipping insertion.")
    
    

def main():
    parser = argparse.ArgumentParser(description="Insert patient pipeline results into the database")
    parser.add_argument("--patient_id", required=True, help="Patient ID for which to insert pipeline results.")
    parser.add_argument("--results_dir", required=True, help="Directory containing the pipeline results files.")
    parser.add_argument("--database_path", required=True, help="Path to the SQLite database file")
    parser.add_argument("--with_update", action='store_true', help="Update existing pipeline results if dataset already exist in the database")
    parser.add_argument("--p_group_file", type=str, required=True, help="Path to P_GROUP_FILE for HLA allele mapping")
    parser.add_argument("--purple_results", type=str, required=True, help="Path to the PURPLE results file, containing ploidy and purity.")
    args = parser.parse_args()
    
    subject_id, glass_nr = args.patient_id.split('-')[0], args.patient_id.split('-')[1]
    
    con = sqlite3.connect(args.database_path)
    cur = con.cursor()
    
    try:
        print(f"Getting current datasets from the database to check for existing records before inserting results for patient {args.patient_id}")
        current_datasets = set((row[0]) for row in cur.execute('SELECT dataset_id FROM allele').fetchall())
        current_datasets = current_datasets.union(set((row[0]) for row in cur.execute('SELECT dataset_id FROM gene_signature_value').fetchall()))
        current_datasets = current_datasets.union(set((row[0]) for row in cur.execute('SELECT dataset_id FROM variant').fetchall()))
        
        print(f"Finding dataset ID for subject_id '{subject_id}' and glass_nr '{glass_nr}'")
        dataset_id = find_dataset_id(con, subject_id, glass_nr)
        print(f"Found dataset ID '{dataset_id}' for subject_id '{subject_id}' and glass_nr '{glass_nr}'")
        
        insert_allele_typing_results(con, args.results_dir, current_datasets, dataset_id, args.p_group_file, with_update=args.with_update)
        insert_allele_loh_results(con, args.results_dir, current_datasets, dataset_id, with_update=args.with_update)
        insert_variant_and_peptide_results(con, args.results_dir, current_datasets, dataset_id, with_update=args.with_update)
        insert_peptide_ms_results(con, args.results_dir, current_datasets, dataset_id, with_update=args.with_update)
        insert_peptide_ms_prediction_results(con, args.results_dir, current_datasets, dataset_id, with_update=args.with_update)
        insert_gene_signature_values(con, args.results_dir, current_datasets, dataset_id, with_update=args.with_update)
        insert_ploidy_and_purity(con, args.purple_results, current_datasets, dataset_id, with_update=args.with_update)
        con.commit()
    except Exception as e:
        print(f"An error occurred while inserting results for patient {args.patient_id}: {e}. Rolling back all changes.")
        con.rollback()
        raise e
    finally:
        con.close()
    
if __name__ == "__main__":
    main()