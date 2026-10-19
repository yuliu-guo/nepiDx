module_path = 'workflow/07_lohhla'
import sys
import pandas as pd

rule prepare_for_lohhla:
    input:
        allele_path = config["results_folder"] + "/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv",
        purity_ploidy_path = config["results_folder"] + "/{patient_id}/07_purple/{patient_id}_tumor_DNA.purple.purity.tsv"
    output:
        allele_path_lohhla = config["results_folder"] + "/{patient_id}/07_lohhla/{patient_id}_tumor_DNA_allele_lohhla.tsv",
        purity_ploidy_lohhla = config["results_folder"] + "/{patient_id}/07_lohhla/{patient_id}_tumor_DNA.purple.purity_lohhla.tsv"
    params:
        prepare_allele_script = module_path + "/scripts/prepare_allele_for_lohhla.py",
        prepare_purity_ploidy_script = module_path + "/scripts/prepare_purity_ploidy_for_lohhla.py",
        sample_name = "{patient_id}",
        p_group_file = config['p_group_file'],
    log:
        config["logs_folder"] + "/run_lohhla/prepare_lohhla_{patient_id}.log"
    shell:
        """
        set -euo pipefail

        mkdir -p $(dirname {log})

        {{
            mkdir -p $(dirname {output.allele_path_lohhla})

            python {params.prepare_allele_script} \
                --p-group-file {params.p_group_file} \
                --input {input.allele_path} \
                --output {output.allele_path_lohhla}

            python {params.prepare_purity_ploidy_script} \
                --input {input.purity_ploidy_path} \
                --output {output.purity_ploidy_lohhla} \
                --sample-name {wildcards.patient_id}
        }} > {log} 2>&1
        """
        
rule run_lohhla:
    input: 
        allele_path_lohhla = rules.prepare_for_lohhla.output.allele_path_lohhla,
        purity_ploidy_lohhla = rules.prepare_for_lohhla.output.purity_ploidy_lohhla,
    output:
        flag = touch(config['results_folder'] + '/flags/07_lohhla/{patient_id}_lohhla.done'),
    log:
        config['logs_folder'] +'/run_lohhla/run_lohhla_{patient_id}.log'
    benchmark:
        'benchmarks/run_lohhla/{patient_id}.log'
    params:
        tmp = get_scratch_directory(['run_lohhla']),
        lohhla_folder = '/ref/lohhla',
        bam_folder_path = (config['scratch'] + '/bam/{patient_id}/'),
        normal_DNA_bam_path = (config['scratch'] + '/bam/{patient_id}/{patient_id}_normal_DNA.bam'),
        #lohhla_sif = (config['singularity_env_folder'] + '/lohhla.sif'),
        lohhla_sif = '/ref/lohhla/lohhla.sif',
        output_dir = config['results_folder'] + '/{patient_id}/07_lohhla',
        HLAexonLoc = '/ref/lohhla/hla.dat',
        novoDir = '/pio/novocraft-V3.09.04',
        minCoverageFilter = 10,
        HLAfastaLoc = '/ref/lohhla/abc_complete.fasta',
        bind_string = lambda wildcards, input, output: get_bind_string(
            input,
            output,
            [],
            get_scratch_directory(['run_lohhla']),
            singularity_internal_dir='/tmp',
            add_singularity_dir=True),
        
    shell:
        """
        (
        mkdir -p {params.output_dir} &&
        mkdir -p {params.tmp} &&

        singularity exec --cleanenv --bind /usr/lib/locale {params.bind_string} \
            --bind {params.lohhla_sif} \
            --bind {params.tmp} \
            --bind {params.bam_folder_path} \
            --bind {params.output_dir} \
            --bind {params.HLAexonLoc}\
            --bind {params.HLAfastaLoc} \
            --bind {params.lohhla_folder} \
            {params.lohhla_sif} \
            Rscript /app/lohhla.R \
            -i {wildcards.patient_id} \
            -o {params.output_dir} \
            -n {params.normal_DNA_bam_path} \
            -b {params.bam_folder_path} \
            -a {input.allele_path_lohhla} \
            -l {params.HLAfastaLoc} \
            -d {params.novoDir} \
            -e {params.HLAexonLoc} \
            -c {input.purity_ploidy_lohhla} \
            -p TRUE \
            -m {params.minCoverageFilter} \
            --fishingStep FALSE \
            --cleanUp FALSE \
        ) > {log} 2>&1
        """
