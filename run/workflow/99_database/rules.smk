module_path = 'workflow/99_database'

rule insert_into_database:
    input:
        hla_allele_path = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
        variants_summary = expand(config['results_folder'] + '/{{patient_id}}/00_variants_summary/{{patient_id}}_{variant}_summary.csv', variant=['SNV', 'genefusion', 'sv']),
        prediction_results = expand(config['results_folder'] + '/{{patient_id}}/03_pMHC_prediction/{variant}_{{patient_id}}_{tool}_with_DAI.csv', tool=['netMHCpan', 'MHCflurry', 'PRIME2', 'netMHCstabpan'], variant=['snv', 'fusion', 'sv']),
        gene_signatures = config['results_folder'] + '/{patient_id}/05_TME/{patient_id}_gene_signatures.csv',
        similarity = expand(config['results_folder'] + '/{{patient_id}}/02_mutated_peptide/{{patient_id}}_{variant}_similarity.csv', variant=['snv', 'fusion', 'sv']),
        immunogenicity = expand(config['results_folder'] + '/{{patient_id}}/08_immunogenicity/{variant}_{{patient_id}}_{tool}.csv', tool=['DeepImmuno', 'ICERFIRE', 'IEDB_immunogenicity'], variant=['snv', 'fusion', 'sv']),
        purple_flag = config['results_folder'] + '/flags/07_purple/{patient_id}_purple.done',
    output:
        flag = touch(config['results_folder'] + '/flags/99_database/{patient_id}.done'),
    threads: 1
    log: config['logs_folder'] +'/99_database/{patient_id}_insert_into_database.log'
    benchmark: 'benchmarks/99_database/{patient_id}_insert_into_database.log'  
    params:
        tmp = get_scratch_directory(['database_insert']),
        data_insert_script = module_path + '/scripts/insert_pipeline_results.py',
        restore_views_script = module_path + '/scripts/create_or_restore_views.py',
        python_for_db_sif = (config['singularity_env_folder'] + '/python_for_db_env.sif'),
        p_group_file = config['p_group_file'],
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [module_path + '/scripts', config['p_group_file']], get_scratch_directory(['database_insert'])),
        results_dir = (config['results_folder'] +'/{patient_id}'),
        database_path = config['patient_database_path'],
        purple_results = config['results_folder'] + '/{patient_id}/07_purple/{patient_id}_tumor_DNA.purple.purity.tsv',
    resources:
        database=1
    shell:
        """
        (
        mkdir -p $(dirname {output.flag}) &&
        mkdir -p {params.tmp} &&
        
        python {params.data_insert_script} \
          --patient_id {wildcards.patient_id} \
          --results_dir {params.results_dir} \
          --database_path {params.database_path} \
          --p_group_file {params.p_group_file} \
          --purple_results {params.purple_results} \
          --with_update &&
          
        python {params.restore_views_script} \
          --database_path {params.database_path}

        ) > {log} 2>&1
        """