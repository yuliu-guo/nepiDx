module_path = 'workflow/00_pre_process'

rule convert_cram_to_bam_tumor_RNA:
    input:
        cram = lambda wc: check_file_exist(wc, "RNA_cram_path"),
    output:
        bam = (config['scratch'] + '/bam/{patient_id}_tumor_RNA.bam'),
        bai = (config['scratch'] + '/bam/{patient_id}_tumor_RNA.bam.bai'),
    threads: 10
    log: config['logs_folder'] +'/convert_cram_to_bam/{patient_id}_tumor_RNA.log'
    benchmark: 'benchmarks/convert_cram_to_bam/{patient_id}_tumor_RNA.log'
    params:
        fasta_ref = config['reference']['GRCh38_fasta'],
    shell:
        """
        (
        module purge
        module load samtools/1.10

        samtools view -b -T {params.fasta_ref} {input.cram} -o {output.bam}
        samtools index -b {output.bam} {output.bai}
        ) > {log} 2>&1


        """
rule convert_cram_to_bam_tumor_DNA:
    input:
        cram = lambda wc: check_file_exist(wc, "Paired_tumor_DNA_cram_path"),
    output:
        bam = (config['scratch'] + '/bam/{patient_id}_tumor_DNA.bam'),
        bai = (config['scratch'] + '/bam/{patient_id}_tumor_DNA.bam.bai'),
    threads: 10
    log: config['logs_folder'] +'/convert_cram_to_bam/{patient_id}_tumor_DNA.log'
    benchmark: 'benchmarks/convert_cram_to_bam/{patient_id}_tumor_DNA.log'
    params:
        fasta_ref = config['reference']['GRCh38_fasta'],
    shell:
        """
        (
        module purge
        module load samtools/1.10

        samtools view -b -T {params.fasta_ref} {input.cram} -o {output.bam}
        samtools index -b {output.bam} {output.bai}
        ) > {log} 2>&1


        """

rule convert_cram_to_bam_normal_DNA:
    input:
        cram = lambda wc: check_file_exist(wc, "Paired_normal_DNA_cram_path"),
    output:
        bam = (config['scratch'] + '/bam/{patient_id}_normal_DNA.bam'),
        bai = (config['scratch'] + '/bam/{patient_id}_normal_DNA.bam.bai'),
    threads: 10
    log: config['logs_folder'] +'/convert_cram_to_bam/{patient_id}_normal_DNA.log'
    benchmark: 'benchmarks/convert_cram_to_bam/{patient_id}_normal_DNA.log'
    params:
        fasta_ref = config['reference']['GRCh38_fasta'],
    shell:
        """
        (
        module purge
        module load samtools/1.10

        samtools view -b -T {params.fasta_ref} {input.cram} -o {output.bam}
        samtools index -b {output.bam} {output.bai}
        ) > {log} 2>&1

        """
