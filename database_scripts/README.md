## Table of Contents
- [Database](#database)

---

## Database

### Directory Structure

```
database_scripts/
├── db_config.json                      # Configuration file for database setup
├── run_db_scripts.py                   # Main script to create the database
├── create_db_and_tables.py             # Script to create the database and the tables
├── insert_datasets_and_signatures.py   # Inserts the list of patients and samples, along with the list of gene signatures
└── get_neoepitope_scores_from_DB.py    # Calculates patient neoepitope scores from the database (check the script for usage instructions)
```

### Creating the Database

Follow these steps to set up the database:

1. Navigate to the `database_scripts` directory:

2. Fill in the correct values, including the path to the database, in `db_config.json`.

3. Run the setup script from the terminal:
   ```bash
   python run_db_scripts.py
   ```

### Inserting Data into the Database

1. Make sure that the path to the database in the Snakemake config is the same as the path to the database in the db_config.json

2. Run the Snakemake pipeline

3. Only patients who have **all required results** from the pipeline will be inserted into the database.


