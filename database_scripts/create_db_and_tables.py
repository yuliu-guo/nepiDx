import sqlite3
import argparse

def create_tables(con):
  cur = con.cursor()

  ### subject table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS subject (
      subject_id    TEXT PRIMARY KEY,
      response      TEXT NOT NULL CHECK (response IN ('CR','PR','PD','SD')),
      sex           TEXT NOT NULL CHECK (sex IN ('male','female')),
      diagnosis     TEXT,
      cohort_name   TEXT NOT NULL,
      os_status     INTEGER,
      os_months     REAL,
      pfs_status     INTEGER,
      pfs_months     REAL,
      created_at    TEXT NOT NULL
    );
  ''')

  ### dataset table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS dataset (
      id            INTEGER PRIMARY KEY,
      subject_id    TEXT NOT NULL,
      glass_nr      TEXT NOT NULL,
      treatment_timeline TEXT NOT NULL CHECK (treatment_timeline IN ('pre','on', 'eot')),
      herv_load     REAL,
      eve_burden    REAL,
      eve_total_tpm  REAL,
      purity        REAL,
      ploidy        REAL,
      TMB           REAL,
      snv_peptide_number     INTEGER,
      snv_variant_number     INTEGER,
      fusion_peptide_number  INTEGER,
      fusion_variant_number  INTEGER,
      sv_peptide_number      INTEGER,
      sv_variant_number      INTEGER,
      created_at    TEXT NOT NULL,

      FOREIGN KEY (subject_id) REFERENCES subject(subject_id),

      UNIQUE (subject_id, glass_nr)
    );
  ''')

  ### variant table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS variant (
      id            INTEGER PRIMARY KEY,
      code          TEXT NOT NULL,
      variant_type  TEXT NOT NULL CHECK (variant_type IN ('SV', 'fusion', 'SNV')),
      dataset_id    INTEGER NOT NULL,
      gene          TEXT,
      created_at    TEXT NOT NULL,

      FOREIGN KEY (dataset_id) REFERENCES dataset(id),

      UNIQUE (dataset_id, code)
    );
  ''')

  ### variant_info table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS variant_info (
      id                INTEGER PRIMARY KEY,
      variant_id        TEXT NOT NULL,
      attribute_name    TEXT NOT NULL,
      attribute_value   TEXT,
      created_at    TEXT NOT NULL,

      FOREIGN KEY (variant_id) REFERENCES variant(id)
    );
  ''')
  
  ### peptide table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS peptide (
      id            INTEGER PRIMARY KEY,
      variant_id    INTEGER NOT NULL,
      sequence      TEXT NOT NULL,
      wildtype_sequence TEXT NOT NULL,
      hamming_similarity  REAL,
      blosum62_2_similarity  REAL,
      created_at    TEXT NOT NULL,

      FOREIGN KEY (variant_id) REFERENCES variant(id)
    );
  ''')
  
  ### peptide_ms table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS peptide_ms (
      id                     INTEGER PRIMARY KEY,
      dataset_id             INTEGER NOT NULL,
      code                   TEXT NOT NULL,
      sequence               TEXT NOT NULL,
      modified_sequence      TEXT NOT NULL,
      n_runs                 INTEGER NOT NULL,
      n_precursors           INTEGER NOT NULL,
      protein_ids            TEXT NOT NULL,
      protein_names          TEXT NOT NULL,
      genes                  TEXT NOT NULL,
      precursor_norm_max     REAL NOT NULL,
      precursor_norm_mean    REAL NOT NULL,
      precursor_norm_median  REAL NOT NULL,
      n_proteins_matched     INTEGER NOT NULL,
      peptide_type           TEXT NOT NULL CHECK (peptide_type IN ('target', 'mutant')),
      created_at             TEXT NOT NULL,

      FOREIGN KEY (dataset_id) REFERENCES dataset(id),
      UNIQUE (dataset_id, code)
    );
  ''')
  

  ### gene_signature_type table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS gene_signature_type (
      name          TEXT PRIMARY KEY,
      gene_list     TEXT NOT NULL,
      method        TEXT NOT NULL
    );
  ''')

  ### gene_signature_value table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS gene_signature_value (
      id                  INTEGER PRIMARY KEY,
      dataset_id          INTEGER NOT NULL,
      gene_signature_type TEXT NOT NULL,
      value               REAL NOT NULL,
      created_at          TEXT NOT NULL,

      FOREIGN KEY (gene_signature_type) REFERENCES gene_signature_type(name),
      FOREIGN KEY (dataset_id) REFERENCES dataset(id)
    );
  ''')

  ### allele table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS allele (
      id            INTEGER PRIMARY KEY,
      dataset_id    INTEGER NOT NULL,
      allele_full_name        TEXT NOT NULL,
      allele_name   TEXT NOT NULL,
      quality       INTEGER NOT NULL,
      created_at    TEXT NOT NULL,

      FOREIGN KEY (dataset_id) REFERENCES dataset(id)
    );
  ''')
  
  ### allele_loh table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS allele_loh (
      id            INTEGER PRIMARY KEY,
      dataset_id    INTEGER NOT NULL,
      HLA           TEXT NOT NULL,
      allele_name_1 TEXT NOT NULL,
      allele_name_2 TEXT NOT NULL,
      copyratio     TEXT NOT NULL,
      kept_HLA      TEXT NOT NULL,
      loss_HLA      TEXT NOT NULL,
      allele_freq_1 REAL NOT NULL,
      allele_freq_2 REAL NOT NULL,
      het_num       INTEGER NOT NULL,
      is_LOH        TEXT NOT NULL CHECK (is_LOH IN ('yes','no')),
      created_at    TEXT NOT NULL,

      FOREIGN KEY (dataset_id) REFERENCES dataset(id)
    );
  ''')

  ### peptide2prediction table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS peptide2prediction (
      id                 INTEGER PRIMARY KEY,
      peptide_id         INTEGER NOT NULL,
      allele_id          INTEGER NOT NULL,
      prediction_tool    TEXT NOT NULL,
      prediction_type    TEXT NOT NULL,
      value              REAL,
      created_at         TEXT NOT NULL,

      FOREIGN KEY (peptide_id)        REFERENCES peptide(id),
      FOREIGN KEY (allele_id)         REFERENCES allele(id),

      UNIQUE (peptide_id, allele_id, prediction_tool, prediction_type)
    );
  ''') 
  
  ### peptide_ms2prediction table
  cur.execute(''' 
    CREATE TABLE IF NOT EXISTS peptide_ms2prediction (
      id                 INTEGER PRIMARY KEY,
      peptide_id         TEXT NOT NULL,
      nmer_sequence      TEXT NOT NULL,
      allele_id          INTEGER NOT NULL,
      prediction_tool    TEXT NOT NULL,
      prediction_type    TEXT NOT NULL,
      value              REAL,
      created_at         TEXT NOT NULL,

      FOREIGN KEY (peptide_id)        REFERENCES peptide_ms(id),
      FOREIGN KEY (allele_id)         REFERENCES allele(id),

      UNIQUE (peptide_id, nmer_sequence, allele_id, prediction_tool, prediction_type)
    );
  ''')  
 

def main():
  parser = argparse.ArgumentParser(description="Insert datasets and gene signatures into the database")
  parser.add_argument("--database_path", required=True, help="Path to the SQLite database file")

  args = parser.parse_args()

  con = sqlite3.connect(args.database_path)

  try:
    create_tables(con)
    con.commit()
    print("Database and tables created.")

  except Exception as e:
    print(f"An error occurred while creating tables: {e}. Rolling back all changes.")
    con.rollback()

  finally:
    con.close()

if __name__ == "__main__":
    main()