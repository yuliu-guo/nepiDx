module_path = 'workflow/03_pMHC_prediction'

rule run_prediction_netMHCpan:
    input: 
        hla_allele_path = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
        mutated_peptides_fasta = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_{variant_type}_nmer.fasta'),
    output:
        netMHCpan_output = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_netMHCpan.csv'),
    log: config['logs_folder'] +'/run_prediction_netMHCpan/{patient_id}_{variant_type}.log'
    benchmark: 'benchmarks/run_prediction_netMHCpan/{patient_id}_{variant_type}.log'  
    params:
        tmp = get_scratch_directory(['netMHCpan_prediction']),
        prediction_netMHCpan = module_path + '/scripts/prediction_netMHCpan.py',
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        p_group_file = config['p_group_file'],
    shell:
        """
        (
        mkdir -p $(dirname {output.netMHCpan_output}) &&
        
        python {params.prediction_netMHCpan} \
            -p {wildcards.patient_id} \
            --hla_allele_path {input.hla_allele_path} \
            --hla_quality_threshold {params.hla_quality_threshold} \
            --mutated_peptides_fasta {input.mutated_peptides_fasta} \
            -o {output.netMHCpan_output} \
            --p_group_file {params.p_group_file}
        ) > {log} 2>&1


        """
    
rule run_prediction_MHCflurry:
    input:
        hla_allele_path = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
        mutated_peptides_fasta = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_{variant_type}_nmer.fasta'),
    output:
        MHCflurry_output = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_MHCflurry.csv'),
    log: config['logs_folder'] +'/run_prediction_MHCflurry/{patient_id}_{variant_type}.log'
    benchmark: 'benchmarks/run_prediction_MHCflurry/{patient_id}_{variant_type}.log'
    params:
        tmp = get_scratch_directory(['MHCflurry_prediction']),
        prediction_mhcflurry = module_path + '/scripts/prediction_MHCflurry.py',
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        p_group_file = config['p_group_file'],
        MHCflurry_sif = (config['singularity_env_folder'] + '/MHCflurry_env.sif'),
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [module_path + '/scripts', config['p_group_file']], get_scratch_directory(['MHCflurry_prediction']), singularity_internal_dir='/tmp', add_singularity_dir=True)
    shell:
        """
        (
        mkdir -p {params.tmp} &&
        
        singularity exec {params.bind_string} -e -H {params.tmp} \
            {params.MHCflurry_sif} \
            python3 {params.prediction_mhcflurry} \
                -p {wildcards.patient_id} \
                --hla_allele_path {input.hla_allele_path} \
                --hla_quality_threshold {params.hla_quality_threshold} \
                --mutated_peptides_fasta {input.mutated_peptides_fasta} \
                -o {output.MHCflurry_output} \
                --scratch_dir {params.tmp} \
                --p_group_file {params.p_group_file}
        ) > {log} 2>&1
        """

rule run_prediction_PRIME2:
    input:
        hla_allele_path = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
        mutated_peptides_fasta = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_{variant_type}_nmer.fasta'),
    output:
        PRIME2_output = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_PRIME2.csv'),
    log: config['logs_folder'] +'/run_prediction_PRIME2/{patient_id}_{variant_type}.log'
    benchmark: 'benchmarks/run_prediction_PRIME2/{patient_id}_{variant_type}.log'
    params:
        tmp = get_scratch_directory(['PRIME2_prediction']),
        prediction_prime2 = module_path + '/scripts/prediction_PRIME2.py',
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        p_group_file = config['p_group_file'],
        PRIME2_sif = (config['singularity_env_folder'] + '/PRIME2_env.sif'),
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [module_path + '/scripts', config['p_group_file']], get_scratch_directory(['PRIME2_prediction']), singularity_internal_dir='/PRIME-2.1/temp', add_singularity_dir=True)
    shell:
        """
        (
        mkdir -p {params.tmp} &&
        
        singularity exec {params.bind_string} -e -H {params.tmp} \
            {params.PRIME2_sif} \
            python3 {params.prediction_prime2} \
                -p {wildcards.patient_id} \
                --hla_allele_path {input.hla_allele_path} \
                --hla_quality_threshold {params.hla_quality_threshold} \
                --mutated_peptides_fasta {input.mutated_peptides_fasta} \
                -o {output.PRIME2_output} \
                --scratch_dir {params.tmp} \
                --p_group_file {params.p_group_file}
                
        ) > {log} 2>&1
        """

rule run_prediction_netMHCstabpan:
    input:
        hla_allele_path = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
        mutated_peptides_fasta = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_{variant_type}_nmer.fasta'),
    output:
        netMHCstabpan_output = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_netMHCstabpan.csv'),
    log: config['logs_folder'] +'/prediction_netMHCstabpan/{patient_id}_{variant_type}.log'
    benchmark: 'benchmarks/run_prediction_netMHCstabpan/{patient_id}_{variant_type}.log'
    params:
        tmp = get_scratch_directory(['netMHCstabpan_prediction']),
        prediction_netMHCstabpan = module_path + '/scripts/prediction_netMHCstabpan.py',
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        p_group_file = config['p_group_file'],
        netMHCstabpan_sif = (config['singularity_env_folder'] + '/netMHCstabpan_env.sif'),
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [module_path + '/scripts', config['p_group_file']], get_scratch_directory(['netMHCstabpan_prediction']), singularity_internal_dir='/tmp', add_singularity_dir=True)
    shell:
        """
        (
        mkdir -p {params.tmp} &&
        
        singularity exec {params.bind_string} -e -H {params.tmp} \
            {params.netMHCstabpan_sif} \
            python3 {params.prediction_netMHCstabpan} \
                -p {wildcards.patient_id} \
                --hla_allele_path {input.hla_allele_path} \
                --hla_quality_threshold {params.hla_quality_threshold} \
                --mutated_peptides_fasta {input.mutated_peptides_fasta} \
                -o {output.netMHCstabpan_output} \
                --scratch_dir {params.tmp} -ia \
                --p_group_file {params.p_group_file}
                
        ) > {log} 2>&1
        """
        
rule run_DAI_netMHCpan:
    input: 
        hla_allele_path = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
        similarity_file = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_{variant_type}_similarity.csv'),
        mutant_results = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_netMHCpan.csv'),
    output:
        DAI_output =  (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_netMHCpan_with_DAI.csv'),
    threads:1
    log: config['logs_folder'] +'/run_prediction_PRIME2/{patient_id}_{variant_type}.log'
    benchmark: 'benchmarks/run_prediction_PRIME2/{patient_id}_{variant_type}.log'  
    params:
        tmp = get_scratch_directory(['netMHCpan_prediction_DAI']),
        DAI_netMHCpan = module_path + '/scripts/DAI_netMHCpan.py',
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        p_group_file = config['p_group_file'],
    shell:
        """
        (
        mkdir -p $(dirname {output.DAI_output}) &&
        mkdir -p {params.tmp} &&

        python {params.DAI_netMHCpan} \
            -p {wildcards.patient_id} \
            --hla_allele_path {input.hla_allele_path} \
            --hla_quality_threshold {params.hla_quality_threshold} \
            --similarity_file {input.similarity_file} \
            --mut_binding_prediction {input.mutant_results} \
            -o {output.DAI_output} \
            --p_group_file {params.p_group_file} \
            --scratch_dir {params.tmp} \
        ) > {log} 2>&1
        """
        
rule run_DAI_MHCflurry:
    input: 
        hla_allele_path = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
        similarity_file = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_{variant_type}_similarity.csv'),
        mutant_results = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_MHCflurry.csv'),
    output:
        DAI_output =  (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_MHCflurry_with_DAI.csv'),
    threads:1
    log: config['logs_folder'] +'/run_DAI_MHCflurry/{patient_id}_{variant_type}.log'
    benchmark: 'benchmarks/run_DAI_MHCflurry/{patient_id}_{variant_type}.log'  
    params:
        tmp = get_scratch_directory(['MHCflurry_prediction_DAI']),
        DAI_MHCflurry = module_path + '/scripts/DAI_MHCflurry.py',
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        p_group_file = config['p_group_file'],
        MHCflurry_sif = (config['singularity_env_folder'] + '/MHCflurry_env.sif'),
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [module_path + '/scripts', config['p_group_file']], get_scratch_directory(['MHCflurry_prediction_DAI']), singularity_internal_dir='/tmp', add_singularity_dir=True)
    shell:
        """
        (
        mkdir -p $(dirname {output.DAI_output}) &&
        mkdir -p {params.tmp} &&
        
        singularity exec {params.bind_string} -e -H {params.tmp} \
            {params.MHCflurry_sif} \
            python3 {params.DAI_MHCflurry} \
                -p {wildcards.patient_id} \
                --hla_allele_path {input.hla_allele_path} \
                --hla_quality_threshold {params.hla_quality_threshold} \
                --similarity_file {input.similarity_file} \
                --mut_binding_prediction {input.mutant_results} \
                -o {output.DAI_output} \
                --p_group_file {params.p_group_file} \
                --scratch_dir {params.tmp} \
        ) > {log} 2>&1
        """
        
rule run_DAI_netMHCstabpan:
    input: 
        hla_allele_path = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
        similarity_file = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_{variant_type}_similarity.csv'),
        mutant_results = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_netMHCstabpan.csv'),
    output:
        DAI_output =  (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_netMHCstabpan_with_DAI.csv'),
    threads:1
    log: config['logs_folder'] +'/run_DAI_netMHCstabpan/{patient_id}_{variant_type}.log'
    benchmark: 'benchmarks/run_DAI_netMHCstabpan/{patient_id}_{variant_type}.log'  
    params:
        tmp = get_scratch_directory(['netMHCstabpan_prediction_DAI']),
        DAI_netMHCstabpan = module_path + '/scripts/DAI_netMHCstabpan.py',
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        p_group_file = config['p_group_file'],
        netMHCstabpan_sif = (config['singularity_env_folder'] + '/netMHCstabpan_env.sif'),
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [module_path + '/scripts', config['p_group_file']], get_scratch_directory(['netMHCstabpan_prediction_DAI']), singularity_internal_dir='/tmp', add_singularity_dir=True)
    shell:
        """
        (
        mkdir -p $(dirname {output.DAI_output}) &&
        mkdir -p {params.tmp} &&
        
        singularity exec {params.bind_string} -e -H {params.tmp} \
            {params.netMHCstabpan_sif} \
            python3 {params.DAI_netMHCstabpan} \
                -p {wildcards.patient_id} \
                --hla_allele_path {input.hla_allele_path} \
                --hla_quality_threshold {params.hla_quality_threshold} \
                --similarity_file {input.similarity_file} \
                --mut_binding_prediction {input.mutant_results} \
                -o {output.DAI_output} -ia \
                --p_group_file {params.p_group_file} \
                --scratch_dir {params.tmp} \
        ) > {log} 2>&1
        """
        
rule run_DAI_PRIME2:
    input: 
        hla_allele_path = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
        similarity_file = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_{variant_type}_similarity.csv'),
        mutant_results = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_PRIME2.csv'),
    output:
        DAI_output =  (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/{variant_type}_{patient_id}_PRIME2_with_DAI.csv'),
    threads:1
    log: config['logs_folder'] +'/run_DAI_PRIME2/{patient_id}_{variant_type}.log'
    benchmark: 'benchmarks/run_DAI_PRIME2/{patient_id}_{variant_type}.log'  
    params:
        tmp = get_scratch_directory(['PRIME2_prediction_DAI']),
        DAI_PRIME2 = module_path + '/scripts/DAI_PRIME2.py',
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        p_group_file = config['p_group_file'],
        PRIME2_sif = (config['singularity_env_folder'] + '/PRIME2_env.sif'),
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [module_path + '/scripts', config['p_group_file']], get_scratch_directory(['PRIME2_prediction_DAI']), singularity_internal_dir='/PRIME-2.1/temp', add_singularity_dir=True)
    shell:
        """
        (
        mkdir -p $(dirname {output.DAI_output}) &&
        mkdir -p {params.tmp} &&
        
        singularity exec {params.bind_string} -e -H {params.tmp} \
            {params.PRIME2_sif} \
            python3 {params.DAI_PRIME2} \
                -p {wildcards.patient_id} \
                --hla_allele_path {input.hla_allele_path} \
                --hla_quality_threshold {params.hla_quality_threshold} \
                --similarity_file {input.similarity_file} \
                --mut_binding_prediction {input.mutant_results} \
                -o {output.DAI_output} \
                --p_group_file {params.p_group_file} \
                --scratch_dir {params.tmp} \
        ) > {log} 2>&1
        """