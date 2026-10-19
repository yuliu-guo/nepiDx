module_path = 'workflow/00_variants_summary'

rule run_SNV_summary:
    input:
        vcf = lambda wc: df.loc[wc.patient_id, 'snv_somatic_annot_vcf_path'],
        cram = lambda wc: check_file_exist(wc, "RNA_cram_path"),
        transcript = lambda wc: check_file_exist(wc, "transcriptome_path"),
        STARcounts = lambda wc: check_file_exist(wc, "STARcounts_RefSeq_path"),
    output:
        snv_summary = (config['results_folder'] + '/{patient_id}/00_variants_summary/{patient_id}_SNV_summary.csv'),
    threads: 20
    resources:
        mem_mb=8000,
        mem_gb = 40,
        walltime_h = 6,
    log: config['logs_folder'] +'/run_summary/{patient_id}.log'
    benchmark: 'benchmarks/run_SNV_summary/{patient_id}.log'  
    params:
        load_file =  module_path + '/scripts/extract_snv_summary.py',
        reference_fasta = config['reference']['GRCh38_fasta'],
    shell:
        """
        (
        module load Conda-BioModules/1.1.0

        mkdir -p $(dirname {output.snv_summary}) &&

        python {params.load_file} \
            --vcf_path {input.vcf} \
            --reference_fasta {params.reference_fasta} \
            --bam_path {input.cram} \
            --STARcounts_path {input.STARcounts} \
            --transcript_path {input.transcript} \
            -q 25 \
            --patient_id {wildcards.patient_id} \
            --output_path {output.snv_summary}
        ) > {log} 2>&1

        """


rule run_genefusion_summary:
    input:
        arriba = lambda wc: check_file_exist(wc, "gene_fusion_arriba"),
        starfusion = lambda wc: check_file_exist(wc, "gene_fusion_star"),
        fusionmap = lambda wc: check_file_exist(wc, "gene_fusion_FusionMap"),
    output:
        fusion_summary = (config['results_folder'] + '/{patient_id}/00_variants_summary/{patient_id}_genefusion_summary.csv'),
    threads: 20
    log: config['logs_folder'] +'/run_genefusion_summary/{patient_id}.log'
    benchmark: 'benchmarks/run_genefusion_summary/{patient_id}.log'  
    resources:
        mem_mb=8000,
        mem_gb=40,
        walltime_h = 6,
    params:
        generate_genefusion_summary =  module_path + '/scripts/generate_genefusion_summary.py',
    shell:
        """
        (
        module load Conda-BioModules/1.1.0
        
        python {params.generate_genefusion_summary} \
            --arriba_file {input.arriba} \
            --starfusion_file {input.starfusion} \
            --fusionmap_file {input.fusionmap} \
            -o {output.fusion_summary}
        ) > {log} 2>&1

        """



rule run_sv_summary:
    input:
        vcf = lambda wc: check_file_exist(wc, "sv_somatic_raw_manta"),
    output:
        out = (config['results_folder'] + '/{patient_id}/00_variants_summary/{patient_id}_sv_summary.csv'),
    threads: 20
    resources:
        mem_mb=8000,
        mem_gb=40,
        walltime_h = 6,
    log: config['logs_folder'] +'/run_sv_summary/{patient_id}.log'
    benchmark: 'benchmarks/run_sv_summary/{patient_id}.log'  
    params:
        script =  module_path + '/scripts/flatten_neosv_allSVTYPE.py',
        kmer = 9,
    shell:
        """
        (
        module load Conda-BioModules/1.1.0

        mkdir -p $(dirname {output.out})

        python {params.script} --vcf {input} --out {output.out} --k {params.kmer}
        
        ) > {log} 2>&1

        """

        

