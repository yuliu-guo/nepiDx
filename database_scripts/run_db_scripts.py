import json
import subprocess
import os

def main():
  with open('db_config.json', 'r') as config_file:
    config_data = json.load(config_file)
  
  print("Running create_db_and_tables.py script.")
  db_directory = os.path.dirname(f"{config_data['database_path']}")
  db_create_bind = f"--bind {db_directory}"
  db_create_command = (
        f"singularity exec {db_create_bind} {config_data['singularity_path']} python3 create_db_and_tables.py --database_path {config_data['database_path']}"
        )
  db_create_result = subprocess.run(db_create_command, check=True, shell=True, stderr=subprocess.STDOUT)
  
  print("Running insert_datasets_and_signatures.py script.")
  data_insert_bind = f"--bind {config_data['database_path']} --bind {config_data['dataset_metadata_filepath']} --bind {config_data['gene_signatures_filepath']} --bind {config_data['herv_load_filepath']}"
  data_insert_command = (
        f"singularity exec {data_insert_bind} --writable-tmpfs {config_data['singularity_path']} python3 insert_datasets_and_signatures.py --database_path {config_data['database_path']} --datasets_file {config_data['dataset_metadata_filepath']} --gene_signatures_file {config_data['gene_signatures_filepath']} --cohort_name {config_data['cohort_name']} --hervload_file {config_data['herv_load_filepath']}"
        )
        
  if config_data['with_update']:
    data_insert_command += " --with_update"
    
  data_insert_result = subprocess.run(data_insert_command, check=True, shell=True, stderr=subprocess.STDOUT)
  
if __name__ == "__main__":
    main()
    