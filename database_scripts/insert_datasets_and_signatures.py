import sqlite3
import pandas as pd
import argparse
import datetime

def get_hervload_from_df(hervload_df, subject_id, glass_nr):
    mask = (hervload_df['subject_id'] == subject_id) & (hervload_df['glass_nr'] == glass_nr)
    filtered_df = hervload_df[mask]
    
    if not filtered_df.empty:
        row = filtered_df.iloc[0]
        return row['herv_load_cpm'], row['eve_burden'], row['eve_total_tpm']
    
    print(f"Subject ID: {subject_id} and glass nr: {glass_nr} not found in the HERV load file. NULL values will be inserted into the database.")
    return None, None, None

def insert_gene_signatures(con, gene_signatures_file, with_update=False):
    """
    Inserts gene signatures data into the database for a given patient.

    Parameters:
        con (sqlite3.Connection): SQLite database connection
        gene_signatures_file (str): Path to the .tsv file containing gene signatures
        with_update (bool): Whether to update existing gene signatures if they already exist in the database (default: False)
    """
    gene_signatures_df = pd.read_csv(gene_signatures_file, sep='\t')
    gene_signatures_df = (
        gene_signatures_df
        .groupby('signatureID')
        .agg(
            gene_list=('Gene', lambda x: ','.join(x)),
            Method=('Method', 'first') 
        )
        .reset_index()
    )
    cur = con.cursor()

    print(f"Starting insert of gene signatures from file '{gene_signatures_file}'")
    existing_signatures = set(row[0] for row in cur.execute('SELECT name FROM gene_signature_type').fetchall())
    for _, row in gene_signatures_df.iterrows():
        if row['signatureID'] not in existing_signatures:
            cur.execute('''
                INSERT INTO gene_signature_type (name, gene_list, method)
                VALUES (?, ?, ?)
            ''', (row['signatureID'], row['gene_list'], row['Method']))
        else:
            if with_update:
                cur.execute('''
                    UPDATE gene_signature_type
                    SET gene_list = ?, method = ?
                    WHERE name = ?
                ''', (row['gene_list'], row['Method'], row['signatureID']))
                print(f"Gene signature '{row['signatureID']}' already exists in the database. Updated gene list to: [{row['gene_list']}] and method to: {row['Method']}.")
            else:
                print(f"Gene signature '{row['signatureID']}' already exists in the database with gene list: [{row['gene_list']}] and method: {row['Method']}. Skipping insertion.")
    print(f"Finished insert of gene signatures from file '{gene_signatures_file}'")

def insert_datasets_and_patients(con, datasets_file, hervload_file, cohort_name, with_update=False):
    """
    Inserts datasets and patient information into the database.

    Parameters:
        con (sqlite3.Connection): SQLite database connection
        datasets_file (str): Path to the .txt file containing datasets and patient information
        hervload_file (str): Path to the .csv file containing HERV load information
        cohort_name (str): The name of the cohort
        with_update (bool): Whether to update existing patient information if they already exist in the database (default: False)
    """
    datasets_df = pd.read_csv(datasets_file, sep='\t', decimal=',', encoding='unicode_escape',
                                usecols=['subject_id', 'Primary_cancer_type',
                                  'os_status', 'os_months', 'pfs_status', 'pfs_months',
                                  'BOR', 'baseline_glass_nr', 'baseline_tmb',
                                  'ontreatment_glass_nr', 'ontreatment_tmb', 'eot_glass'],
                                dtype={'subject_id': str, 'Primary_cancer_type': str,
                                       'os_status': int, 'os_months': float, 
                                       'pfs_status': int, 'pfs_months': float,
                                       'BOR': str, 'baseline_glass_nr': str, 'baseline_tmb': str,
                                       'ontreatment_glass_nr': str, 'ontreatment_tmb': str,
                                       'eot_glass': str})
                                    
    datasets_df['baseline_tmb'] = pd.to_numeric(datasets_df['baseline_tmb'], downcast="float")
    datasets_df['ontreatment_tmb'] = pd.to_numeric(datasets_df['ontreatment_tmb'], downcast="float")
    hervload_df = pd.read_csv(hervload_file, sep=",")
    hervload_df['subject_id'] = hervload_df['patient_id'].str.split("-").str[0]
    hervload_df['glass_nr'] = hervload_df['patient_id'].str.split("-").str[1]
    hervload_df = hervload_df.drop(columns=['Primary_cancer_type', 'patient_id'])
    
    cur = con.cursor()

    ### Insert subjects
    print(f"Starting insert of patients from file '{datasets_file}'")
    current_subjects = set(row[0] for row in cur.execute('SELECT subject_id FROM subject').fetchall())
    datasets_df['sex'] = datasets_df['subject_id'].apply(lambda x: 'male' if x.endswith('m') else 'female')
    for _, row in datasets_df.iterrows():
        if row['subject_id'] not in current_subjects:
            cur.execute('''
                INSERT INTO subject (subject_id, response, sex, diagnosis, cohort_name, os_status, os_months, pfs_status, pfs_months, created_at)
                VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ''', (row['subject_id'], row['BOR'], 
                  row['sex'], row['Primary_cancer_type'], 
                  cohort_name, row['os_status'], row['os_months'],
                  row['pfs_status'], row['pfs_months'],
                  datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        else:
            if with_update:
                cur.execute('''
                    UPDATE subject
                    SET response = ?, sex = ?, diagnosis = ?, cohort_name = ?, 
                        os_status = ?, os_months = ?, pfs_status = ?, pfs_months = ?, 
                        created_at = ?
                    WHERE subject_id = ?
                ''', (row['BOR'], row['sex'], row['Primary_cancer_type'], cohort_name,
                      row['os_status'], row['os_months'],
                      row['pfs_status'], row['pfs_months'],
                     datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), row['subject_id']))
                     
                print(f"Subject '{row['subject_id']}' already exists in the database. Updated response to: {row['BOR']}, sex to: {row['sex']}, diagnosis to: {row['Primary_cancer_type']} and cohort to: {cohort_name}.")
            else:
                print(f"Subject '{row['subject_id']}' already exists in the database with response: {row['BOR']}, sex: {row['sex']}, diagnosis: {row['Primary_cancer_type']} and cohort: {cohort_name}. Skipping insertion.")
    print(f"Finished insert of patients from file '{datasets_file}'")
    
    current_datasets = set((row[0], row[1]) for row in cur.execute('SELECT subject_id, glass_nr FROM dataset').fetchall())
    
    ### Insert baseline datasets
    print(f"Starting insert of datasets from file '{datasets_file}'")
    for _, row in datasets_df.iterrows():
        if row['baseline_glass_nr'] and str(row['baseline_glass_nr']).strip() != '' and str(row['baseline_glass_nr']).lower() != 'nan':
        
            ### Find HERV load for the given dataset
            hervload, eve_burden, eve_tpm = get_hervload_from_df(hervload_df, row['subject_id'], str(row['baseline_glass_nr']))
            
            if (row['subject_id'], str(row['baseline_glass_nr'])) not in current_datasets:
                cur.execute('''
                    INSERT INTO dataset (subject_id, glass_nr, treatment_timeline,
                                         herv_load, eve_burden, eve_total_tpm, TMB, created_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                ''', (row['subject_id'], str(row['baseline_glass_nr']), 'pre', hervload, eve_burden, eve_tpm, row['baseline_tmb'], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            else:
                if with_update:
                    cur.execute('''
                        UPDATE dataset
                        SET treatment_timeline = ?, herv_load = ?, eve_burden = ?,
                            eve_total_tpm = ?, TMB = ?, created_at = ?
                        WHERE subject_id = ? AND glass_nr = ?
                    ''', ('pre', hervload, eve_burden, eve_tpm, row['baseline_tmb'], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), row['subject_id'], str(row['baseline_glass_nr'])))
                    
                    print(f"Dataset for subject '{row['subject_id']}' with glass number '{str(row['baseline_glass_nr']).replace('-', '')}' already exists in the database. Updated treatment timeline to: 'pre', HERV load to: {hervload}, EVE burden to: {eve_burden}, EVE total TPM to: {eve_tpm}, and TMB to: {row['baseline_tmb']}.")
                else:
                    print(f"Dataset for subject '{row['subject_id']}' with glass number '{str(row['baseline_glass_nr']).replace('-', '')}' already exists in the database. Skipping insertion.")

    ### Insert ontreatment datasets
    for _, row in datasets_df.iterrows():
        if row['ontreatment_glass_nr'] and str(row['ontreatment_glass_nr']).strip() != '' and str(row['ontreatment_glass_nr']).lower() != 'nan':
            
            ### Find HERV load for the given dataset
            hervload, eve_burden, eve_tpm = get_hervload_from_df(hervload_df, row['subject_id'], str(row['ontreatment_glass_nr']))
            
            if (row['subject_id'], str(row['ontreatment_glass_nr'])) not in current_datasets:
                cur.execute('''
                    INSERT INTO dataset (subject_id, glass_nr, treatment_timeline,
                                         herv_load, eve_burden, eve_total_tpm, TMB, created_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                ''', (row['subject_id'], str(row['ontreatment_glass_nr']), 'on',  hervload, eve_burden, eve_tpm, row['ontreatment_tmb'], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            else:
                if with_update:
                    cur.execute('''
                        UPDATE dataset
                        SET treatment_timeline = ?, herv_load = ?, eve_burden = ?,
                            eve_total_tpm = ?, TMB = ?, created_at = ?
                        WHERE subject_id = ? AND glass_nr = ?
                    ''', ('on', hervload, eve_burden, eve_tpm, row['ontreatment_tmb'], datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), row['subject_id'], str(row['ontreatment_glass_nr'])))
                    
                    print(f"Dataset for subject '{row['subject_id']}' with glass number '{str(row['ontreatment_glass_nr'])}' already exists in the database. Updated treatment timeline to: 'on', HERV load to: {hervload}, EVE burden to: {eve_burden}, EVE total TPM to: {eve_tpm}, and TMB to: {row['ontreatment_tmb']}.")
                else:
                    print(f"Dataset for subject '{row['subject_id']}' with glass number '{str(row['ontreatment_glass_nr'])}' already exists in the database. Skipping insertion.")

    ### Insert EOT datasets
    for _, row in datasets_df.iterrows():
        if row['eot_glass'] and str(row['eot_glass']).strip() != '' and str(row['eot_glass']).lower() != 'nan':
            
            ### Find HERV load for the given dataset
            hervload, eve_burden, eve_tpm = get_hervload_from_df(hervload_df, row['subject_id'], str(row['eot_glass']).replace('-', ''))
            
            if (row['subject_id'], str(row['eot_glass']).replace('-', '')) not in current_datasets:
                cur.execute('''
                    INSERT INTO dataset (subject_id, glass_nr, treatment_timeline,
                                         herv_load, eve_burden, eve_total_tpm, created_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                ''', (row['subject_id'], str(row['eot_glass']).replace('-', ''), 'eot',  hervload, eve_burden, eve_tpm, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
            else:
                if with_update:
                    cur.execute('''
                        UPDATE dataset
                        SET treatment_timeline = ?, herv_load = ?, eve_burden = ?,
                            eve_total_tpm = ?, created_at = ?
                        WHERE subject_id = ? AND glass_nr = ?
                    ''', ('eot', hervload, eve_burden, eve_tpm, datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"), row['subject_id'], str(row['eot_glass']).replace('-', '')))
                    
                    print(f"Dataset for subject '{row['subject_id']}' with glass number '{str(row['eot_glass']).replace('-', '')}' already exists in the database. Updated treatment timeline to: 'eot', HERV load to: {hervload}, EVE burden to: {eve_burden}, and EVE total TPM to: {eve_tpm}.")
                else:
                    print(f"Dataset for subject '{row['subject_id']}' with glass number '{str(row['eot_glass']).replace('-', '')}' already exists in the database. Skipping insertion.")
    print(f"Finished insert of datasets from file '{datasets_file}'")

def main():
    parser = argparse.ArgumentParser(description="Insert datasets and gene signatures into the database")
    parser.add_argument("--gene_signatures_file", help="Path to the gene signatures Excel file")
    parser.add_argument("--datasets_file", help="Path to the datasets .txt file")
    parser.add_argument("--hervload_file", help="Path to the hervload .csv file")
    parser.add_argument("--cohort_name", help="The name of the cohort")
    parser.add_argument("--database_path", required=True, help="Path to the SQLite database file")
    parser.add_argument("--with_update", action='store_true', help="Update existing gene signatures if they already exist in the database")

    args = parser.parse_args()

    con = sqlite3.connect(args.database_path)
    
    try:
        if args.gene_signatures_file:
            insert_gene_signatures(con, args.gene_signatures_file, with_update=args.with_update)

        if args.datasets_file:
            if args.cohort_name and args.hervload_file:
                insert_datasets_and_patients(con, args.datasets_file, args.hervload_file, args.cohort_name, with_update=args.with_update)
            else:
                raise Exception('Provide a datasets file, a HERV load file, and a cohort name to insert patients and datasets into the database.')
        con.commit()

    except Exception as e:
        print(f"An error occurred when inserting datasets and gene signatures: {e}. Rolling back all changes.")
        con.rollback()
        raise e

    finally:
        con.close()

if __name__ == "__main__":
    main()