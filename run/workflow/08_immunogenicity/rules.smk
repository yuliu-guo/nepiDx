module_path = 'workflow/08_immunogenicity'

rule DeepImmuno_prediction:
    input:
        peptides = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_{variant_type}_nmer.fasta'),
        alleles = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
    output:
        predictions = (config['results_folder'] + '/{patient_id}/08_immunogenicity/{variant_type}_{patient_id}_DeepImmuno.csv'),
    threads:1
    log: config['logs_folder'] +'/prediction_DeepImmuno/{patient_id}_{variant_type}.log'
    benchmark: 'benchmarks/run_prediction_DeepImmuno/{patient_id}_{variant_type}.log'
    params:
        tmp = get_scratch_directory(['DeepImmuno_prediction']),
        deepimmuno_sif = (config['singularity_env_folder'] + '/DeepImmuno_env.sif'),
        prediction_script = module_path + '/scripts/DeepImmuno_predict.py',
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [module_path + '/scripts', config['p_group_file']], get_scratch_directory(['DeepImmuno_prediction'])),
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        p_group_file = config['p_group_file'],
    shell:
        """
        (
        mkdir -p $(dirname {output.predictions}) &&
        mkdir -p {params.tmp} &&

        singularity run {params.bind_string} -e -H {params.tmp} \
        {params.deepimmuno_sif} \
            python3 {params.prediction_script} \
            --input_peptides {input.peptides} \
            --input_alleles {input.alleles} \
            --output_file {output.predictions} \
            --scratch_dir {params.tmp} \
            --patient_id {wildcards.patient_id} \
            --p_group_file {params.p_group_file} \
            --hla_quality_threshold {params.hla_quality_threshold}

        ) > {log} 2>&1
        """

rule IEDB_prediction:
    input:
        peptides = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_{variant_type}_nmer.fasta'),
        alleles = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
    output:
        predictions = (config['results_folder'] + '/{patient_id}/08_immunogenicity/{variant_type}_{patient_id}_IEDB_immunogenicity.csv'),
    threads:1
    log: config['logs_folder'] +'/prediction_IEDB/{patient_id}_{variant_type}.log'
    benchmark: 'benchmarks/run_prediction_IEDB/{patient_id}_{variant_type}.log'
    params:
        tmp = get_scratch_directory(['IEDB_prediction']),
        prediction_script = module_path + '/scripts/IEDB_predict.py',
        alleles_list = module_path + '/scripts/IEDB_allele_list.txt',
        IEDB_sif = (config['singularity_env_folder'] + '/IEDB_immunogenicity_env.sif'),
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [module_path + '/scripts', config['p_group_file']], get_scratch_directory(['IEDB_prediction'])),
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        p_group_file = config['p_group_file'],
    shell:
        """
        (
        mkdir -p $(dirname {output.predictions}) &&
        mkdir -p {params.tmp} &&

        singularity exec {params.bind_string} -e -H {params.tmp} \
        {params.IEDB_sif} \
            python3 {params.prediction_script} \
                --input_peptides {input.peptides} \
                --input_alleles {input.alleles} \
                --output_file {output.predictions} \
                --alleles_list {params.alleles_list} \
                --scratch_dir {params.tmp} \
                --patient_id {wildcards.patient_id} \
                --p_group_file {params.p_group_file} \
                --hla_quality_threshold {params.hla_quality_threshold}


        ) > {log} 2>&1
        """

rule ICERFIRE_prediction:
    input:
        peptides = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_{variant_type}_similarity.csv'),
        alleles = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
    output:
        predictions = (config['results_folder'] + '/{patient_id}/08_immunogenicity/{variant_type}_{patient_id}_ICERFIRE.csv'),
    log: config['logs_folder'] +'/prediction_ICERFIRE/{patient_id}_{variant_type}.log'
    benchmark: 'benchmarks/run_prediction_ICERFIRE/{patient_id}_{variant_type}.log'
    params:
        tmp = get_scratch_directory(['ICERFIRE_prediction']),
        prediction_script = module_path + '/scripts/ICERFIRE_predict.py',
        alleles_list = module_path + '/scripts/ICERFIRE_allele_list.txt',
        ICERFIRE_sif = (config['singularity_env_folder'] + '/ICERFIRE_no_expr_env.sif'),
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [module_path + '/scripts', config['p_group_file']], get_scratch_directory(['ICERFIRE_prediction'])),
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        p_group_file = config['p_group_file'],
    shell:
        """
        (
        mkdir -p $(dirname {output.predictions}) &&
        mkdir -p {params.tmp} &&

        singularity exec --writable-tmpfs -B {params.tmp}:/ICERFIRE/tmp/ \
        -B {params.tmp}:/ICERFIRE/output/ {params.bind_string} -e -H {params.tmp} \
        {params.ICERFIRE_sif} \
            python3 {params.prediction_script} \
                --input_peptides {input.peptides} \
                --input_alleles {input.alleles} \
                --output_file {output.predictions} \
                --alleles_list {params.alleles_list} \
                --scratch_dir {params.tmp} \
                --patient_id {wildcards.patient_id} \
                --p_group_file {params.p_group_file} \
                --hla_quality_threshold {params.hla_quality_threshold}

        ) > {log} 2>&1
        """