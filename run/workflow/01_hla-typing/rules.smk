module_path = 'workflow/01_hla-typing'


rule t1k_tumor_RNA:
    input: 
        bam = rules.convert_cram_to_bam_tumor_RNA.output.bam,
    output:
        allele = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_RNA_allele.tsv'),
    log: config['logs_folder'] +'/hla-typing/{patient_id}_tumor_RNA.log'
    threads: 20
    resources:
        mem_mb=8000,
        mem_gb=40,
        walltime_h = 6,
    benchmark: 'benchmarks/hla_typing_t1k/{patient_id}_tumor_RNA.log'
    params:
        coord = config['reference']['hla_rna_coord'],
        hla_rna_ref = config['reference']['hla_rna_ref'],
        tmp = get_scratch_directory(['hla-typing']),
    shell:
        r"""
        mkdir -p $(dirname {output.allele})
        (module load t1k/1.0.8
        run-t1k -t 35 -c {params.coord} -f {params.hla_rna_ref} -b {input} -o {wildcards.patient_id} --od {params.tmp} --preset hla
            
        cp {params.tmp}/{wildcards.patient_id}_allele.tsv {output.allele}
        ) > {log} 2>&1
           
        """


rule t1k_tumor_DNA:
    input: 
        bam = rules.convert_cram_to_bam_tumor_DNA.output.bam,
    output:
        allele = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_tumor_DNA_allele.tsv'),
    log: config['logs_folder'] +'/hla-typing/{patient_id}_tumor_DNA.log'
    threads: 20
    resources:
        mem_mb=8000,
        mem_gb=40,
        walltime_h = 6,
    benchmark: 'benchmarks/hla_typing_t1k/{patient_id}_tumor_DNA.log'
    params:
        coord = config['reference']['hla_dna_coord'],
        hla_dna_ref = config['reference']['hla_dna_ref'],
        tmp = get_scratch_directory(['hla-typing']),
    shell:
        r"""
        mkdir -p $(dirname {output.allele})
        (module load t1k/1.0.8
        run-t1k -t 35 -c {params.coord} -f {params.hla_dna_ref} -b {input} -o {wildcards.patient_id} --od {params.tmp} --preset hla-wgs
            
        cp {params.tmp}/{wildcards.patient_id}_allele.tsv {output.allele}
        ) > {log} 2>&1
           
        """


rule t1k_normal_DNA:
    input: 
        bam = rules.convert_cram_to_bam_normal_DNA.output.bam,
    output:
        allele = (config['results_folder'] +'/{patient_id}/01_hla-typing/{patient_id}_normal_DNA_allele.tsv'),
    log: config['logs_folder'] +'/hla-typing/{patient_id}_normal_DNA.log'
    threads: 20
    resources:
        mem_mb=8000,
        mem_gb=40,
        walltime_h = 6,
    benchmark: 'benchmarks/hla_typing_t1k/{patient_id}_normal_DNA.log'
    params:
        coord = config['reference']['hla_dna_coord'],
        hla_dna_ref = config['reference']['hla_dna_ref'],
        tmp = get_scratch_directory(['hla-typing']),
    shell:
        r"""
        mkdir -p $(dirname {output.allele})
        (module load t1k/1.0.8
        run-t1k -t 35 -c {params.coord} -f {params.hla_dna_ref} -b {input} -o {wildcards.patient_id} --od {params.tmp} --preset hla-wgs
            
        cp {params.tmp}/{wildcards.patient_id}_allele.tsv {output.allele}
        ) > {log} 2>&1
           
        """