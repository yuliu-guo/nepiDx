module_path = 'workflow/06_isovar/'

rule run_isovar:
    input: 
        #unpack(get_isovar_input)
        vcf = lambda wc: df.loc[wc.patient_id, 'snv_somatic_annot_vcf_path'],
        #bam = rules.convert_cram_to_bam_tumor_RNA.output.bam,
        bam = (config['scratch'] + '/bam/{patient_id}_tumor_RNA.bam'),
        
    output:
        csv = (config['results_folder'] +'/{patient_id}/06_isovar/{patient_id}_isovar.csv'),
    threads: 10
    log: config['logs_folder'] +'/run_isovar/{patient_id}.log'
    benchmark: 'benchmarks/run_isovar/{patient_id}.log'
    params:
        tmp = get_scratch_directory(['isovar']),
        #reference_fasta = config['reference']['GRCh38_fasta'],
        isovar_sif = '/ref/isovar_env.sif',
        bind_string = lambda wildcards, input, output: get_bind_string(input, output, [], get_scratch_directory(['isovar']), singularity_internal_dir='/tmp', add_singularity_dir=True)
    shell:
        """
        (
        mkdir -p {params.tmp} &&

        tmp_vcf="$(mktemp -p {params.tmp} isovar.vcf.XXXXXX)" &&
        trap 'rm -f "$tmp_vcf"' EXIT &&

        gzip -cd "{input.vcf}" | sed '/^##FILTER=<ID=gnomAD_AF_too_high,/d' > "$tmp_vcf" &&

        singularity exec {params.bind_string} -e -H {params.tmp} \
            {params.isovar_sif} isovar --vcf "$tmp_vcf" \
            --genome GRCh38 \
            --bam {input.bam} \
            --output {output.csv}
        ) > {log} 2>&1

        
        """
        