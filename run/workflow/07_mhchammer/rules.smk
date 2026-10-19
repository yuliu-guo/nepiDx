module_path = 'workflow/07_mhchammer'
import sys
import pandas as pd

rule run_amber:
    input:
        normal_wxs_bam = (config['scratch'] + '/bam/{patient_id}_normal_DNA.bam'),
        normal_wxs_bai = (config['scratch'] + '/bam/{patient_id}_normal_DNA.bam.bai'),
        tumor_wxs_bam  = (config['scratch'] + '/bam/{patient_id}_tumor_DNA.bam'),
        tumor_wxs_bai  = (config['scratch'] + '/bam/{patient_id}_tumor_DNA.bam.bai'),
    output:
        flag = touch(config['results_folder'] + '/flags/07_amber/{patient_id}_amber.done'),
    log:
        config['logs_folder'] +'/run_amber/{patient_id}.log'
    benchmark:
        'benchmarks/run_amber/{patient_id}.log'
    params:
        tmp = get_scratch_directory(['run_amber']),
        amber_sif = (config['singularity_env_folder'] + '/amber_4.2.sif'),
        output_dir = config['results_folder'] + '/{patient_id}/07_amber',
        bind_string = lambda wildcards, input, output: get_bind_string(
            input,
            output,
            [],
            get_scratch_directory(['run_amber']),
            singularity_internal_dir='/tmp',
            add_singularity_dir=True
        ),
        reference_genome =  config['reference']['GRCh38_fasta'],
        loci_file = config['purple_ref_dir'] + '/AmberGermlineSites.38.tsv.gz',
        threads = 16,
    shell:
        """
        (
        mkdir -p {params.output_dir} &&
        mkdir -p {params.tmp} &&

        cp {input.normal_wxs_bam} {input.normal_wxs_bai} \
           {input.tumor_wxs_bam} {input.tumor_wxs_bai} \
           {params.tmp} &&

        singularity run --bind /usr/lib/locale {params.bind_string} \
            --bind {params.amber_sif} \
            --bind {params.reference_genome} \
            --bind {params.loci_file} \
            {params.amber_sif} \
            -reference {wildcards.patient_id}_normal_DNA \
            -reference_bam {params.tmp}/{wildcards.patient_id}_tumor_DNA.bam \
            -tumor {wildcards.patient_id}_tumor_DNA \
            -tumor_bam {params.tmp}/{wildcards.patient_id}_tumor_DNA.bam \
            -output_dir {params.output_dir}/ \
            -threads {params.threads} \
            -loci {params.loci_file} \
            -ref_genome {params.reference_genome} \
            -ref_genome_version 38 
        ) > {log} 2>&1
        """

rule run_cobalt:
    input:
        normal_wxs_bam = (config['scratch'] + '/bam/{patient_id}_normal_DNA.bam'),
        normal_wxs_bai = (config['scratch'] + '/bam/{patient_id}_normal_DNA.bam.bai'),
        tumor_wxs_bam  = (config['scratch'] + '/bam/{patient_id}_tumor_DNA.bam'),
        tumor_wxs_bai  = (config['scratch'] + '/bam/{patient_id}_tumor_DNA.bam.bai'),
    output:
        flag = touch(config['results_folder'] + '/flags/07_cobalt/{patient_id}_cobalt.done'),
    log:
        config['logs_folder'] +'/run_cobalt/{patient_id}.log'
    benchmark:
        'benchmarks/run_cobalt/{patient_id}.log'
    params:
        tmp = get_scratch_directory(['run_cobalt']),
        cobalt_sif = (config['singularity_env_folder'] + '/cobalt_2.2.sif'),
        output_dir = config['results_folder'] + '/{patient_id}/07_cobalt',
        bind_string = lambda wildcards, input, output: get_bind_string(
            input,
            output,
            [],
            get_scratch_directory(['run_cobalt']),
            singularity_internal_dir='/tmp',
            add_singularity_dir=True
        ),
        reference_genome = config['reference']['GRCh38_fasta'],
        gc_profile = config['purple_ref_dir'] + '/GC_profile.1000bp.38.cnp',
        threads = 16,
    shell:
        """
        (
        mkdir -p {params.output_dir} &&
        mkdir -p {params.tmp} &&

        cp {input.normal_wxs_bam} {input.normal_wxs_bai} \
           {input.tumor_wxs_bam} {input.tumor_wxs_bai} \
           {params.tmp} &&

        singularity run --bind /usr/lib/locale {params.bind_string} \
            --bind {params.cobalt_sif} \
            --bind {params.reference_genome} \
            --bind {params.gc_profile} \
            {params.cobalt_sif} \
            -reference {wildcards.patient_id}_normal_DNA \
            -reference_bam {params.tmp}/{wildcards.patient_id}_normal_DNA.bam \
            -tumor {wildcards.patient_id}_tumor_DNA \
            -tumor_bam {params.tmp}/{wildcards.patient_id}_tumor_DNA.bam \
            -output_dir {params.output_dir}/ \
            -threads {params.threads} \
            -gc_profile {params.gc_profile} \
            -ref_genome {params.reference_genome} \
            -ref_genome_version 38
        ) > {log} 2>&1
        """

rule run_purple:
    input:
        amber_flag  = config['results_folder'] + '/flags/07_amber/{patient_id}_amber.done',
        cobalt_flag = config['results_folder'] + '/flags/07_cobalt/{patient_id}_cobalt.done',
    output:
        flag = touch(config['results_folder'] + '/flags/07_purple/{patient_id}_purple.done'),
    log:
        config['logs_folder'] +'/run_purple/{patient_id}.log'
    benchmark:
        'benchmarks/run_purple/{patient_id}.log'
    params:
        tmp = get_scratch_directory(['run_purple']),
        purple_sif = (config['singularity_env_folder'] + '/purple_4.3.sif'),
        amber_dir = config['results_folder'] + '/{patient_id}/07_amber',
        cobalt_dir = config['results_folder'] + '/{patient_id}/07_cobalt',
        output_dir = config['results_folder'] + '/{patient_id}/07_purple',
        bind_string = lambda wildcards, input, output: get_bind_string(
            input,
            output,
            [],
            get_scratch_directory(['run_purple']),
            singularity_internal_dir='/tmp',
            add_singularity_dir=True),
        reference_genome = config['reference']['GRCh38_fasta'],
        gc_profile = config['purple_ref_dir'] + '/GC_profile.1000bp.38.cnp',
        ensembl_data_cache = config['purple_ref_dir'] + '/ensembl_data_cache',
        driver_gene_panel = config['purple_ref_dir'] + '/DriverGenePanel.38.tsv',
        somatic_hotspots = config['purple_ref_dir'] + '/KnownHotspots.somatic.38.vcf.gz',
        germline_hotspots = config['purple_ref_dir'] + '/KnownHotspots.germline.38.vcf.gz',
        germline_del_freq = config['purple_ref_dir'] + '/cohort_germline_del_freq.38.csv',
        circos_bin = '/opt/circos/circos-0.69-9/bin/circos',
    shell:
        """
        (
        mkdir -p {params.output_dir} &&
        mkdir -p {params.tmp} &&

        singularity run --bind /usr/lib/locale {params.bind_string} \
            --bind {params.purple_sif} \
            --bind {params.amber_dir} \
            --bind {params.cobalt_dir} \
            --bind {params.reference_genome} \
            --bind {params.reference_genome}.fai \
            --bind {params.gc_profile} \
            --bind {params.ensembl_data_cache} \
            --bind {params.driver_gene_panel} \
            --bind {params.somatic_hotspots} \
            --bind {params.germline_hotspots} \
            --bind {params.germline_del_freq} \
            {params.purple_sif} \
            -reference {wildcards.patient_id}_normal_DNA \
            -tumor {wildcards.patient_id}_tumor_DNA \
            -amber {params.amber_dir} \
            -cobalt {params.cobalt_dir} \
            -ref_genome_version 38 \
            -ref_genome {params.reference_genome} \
            -gc_profile {params.gc_profile} \
            -ensembl_data_dir {params.ensembl_data_cache}/ \
            -driver_gene_panel {params.driver_gene_panel} \
            -somatic_hotspots {params.somatic_hotspots} \
            -germline_hotspots {params.germline_hotspots} \
            -germline_del_freq_file {params.germline_del_freq} \
            -circos {params.circos_bin} \
            -output_dir {params.output_dir}/
        ) > {log} 2>&1
        """


rule run_mhchammer:
    input:
        tumor_rna_bam = (config['scratch'] + '/bam/{patient_id}_RNA.bam'),
        #tumor_rna_bam = rules.convert_cram_to_bam_tumor_RNA.output.bam,
        #normal_wxs_bam = rules.convert_cram_to_bam_normal_DNA.output.bam,
        normal_wxs_bam = (config['scratch'] + '/bam/{patient_id}_normal_DNA.bam'),
        tumor_wxs_bam = (config['scratch'] + '/bam/{patient_id}_tumor_DNA.bam'),
        #tumor_wxs_bam = rules.convert_cram_to_bam_tumor_DNA.output.bam,
        #allele_path = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_RNA_allele.tsv'),
        purple_flag = config['results_folder'] + '/flags/07_purple/{patient_id}_purple.done',
    output:
        flag = touch(config['results_folder'] + '/flags/07_mhchammer/{patient_id}.done'),
    log: config['logs_folder'] +'/mhchammer/run_mhchammer_{patient_id}.log'
    benchmark: 'benchmarks/run_mhchammer/{patient_id}.log'
    params:
        tmp = get_scratch_directory(['run_mhchammer']),
        mhc_hammer_inventory_script = module_path + '/scripts/create_mhchammer_inventory.py',
        mhchammer_sif = (config['singularity_env_folder'] + '/mhchammer_env.sif'),
        output_dir = (config['results_folder'] + '/{patient_id}/07_mhchammer'),
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [], get_scratch_directory(['run_mhchammer']), singularity_internal_dir='/mhc-hammer/mhc_hammer_temp_dir', add_singularity_dir=True),
        nextflow_log_file = 'logs/mhchammer/.nextflow.log',
        hla_quality_threshold = config['params']['hla_quality_threshold'],
        nextflow_config_file = module_path + '/scripts/nextflow.hpc.config',
        purple_results = config['results_folder'] + '/{patient_id}/07_purple/{patient_id}_tumor_DNA.purple.purity.tsv',
    handover: True
    shell:
        """
        (
        mkdir -p {params.output_dir} &&
        mkdir -p {params.tmp} &&
        cp {input.tumor_rna_bam} {input.normal_wxs_bam} {input.tumor_wxs_bam} \
         {input.tumor_rna_bam}.bai {input.normal_wxs_bam}.bai {input.tumor_wxs_bam}.bai \
         {params.purple_results} \
         {params.tmp} &&
        
        singularity exec {params.bind_string} -e -H {params.tmp} {params.mhchammer_sif} \
            python3 {params.mhc_hammer_inventory_script} \
            -p {wildcards.patient_id} \
            --normal_wxs {params.tmp}/$(basename {input.normal_wxs_bam}) \
            --tumor_wxs {params.tmp}/$(basename {input.tumor_wxs_bam}) \
            --tumor_rna {params.tmp}/$(basename {input.tumor_rna_bam}) \
            --purple_result_file {params.tmp}/$(basename {params.purple_results}) \
            --quality_threshold {params.hla_quality_threshold} \
            --scratch_dir {params.tmp} &&
            
        singularity exec --writable-tmpfs {params.bind_string} -e -H {params.tmp} \
        --env NXF_HOME=/nextflow-exec --env NXF_OFFLINE='true' --env NXF_LOG_FILE={params.nextflow_log_file} \
        {params.mhchammer_sif} \
         nextflow run /mhc-hammer/main.nf -profile singularity -c {params.nextflow_config_file} \
          --input {params.tmp}/{wildcards.patient_id}_mhchammer_inventory.csv --hlahd_local_install false \
          --outdir {params.output_dir} --tracedir $(dirname {log}) &&
          
        rm -f {params.tmp}/{wildcards.patient_id}_mhchammer_inventory.csv {params.tmp}/{wildcards.patient_id}_alleles.csv \
         {params.tmp}/$(basename {input.normal_wxs_bam}) {params.tmp}/$(basename {input.tumor_wxs_bam}) \
         {params.tmp}/$(basename {input.tumor_rna_bam}) {params.tmp}/$(basename {input.normal_wxs_bam}).bai \
         {params.tmp}/$(basename {input.tumor_wxs_bam}).bai {params.tmp}/$(basename {input.tumor_rna_bam}).bai
        
        ) > {log} 2>&1
        """

