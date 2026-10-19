module_path = 'workflow/05_TME'

    
rule gene_signature_scores:
    input:
        #gene_counts = lambda wc: df.loc[wc.patient_id, 'STARcounts_RefSeq_path'],
        gene_counts = lambda wc: check_file_exist(wc, "STARcounts_RefSeq_path"),
        #gene_signature_scores_input
        #lambda wc: df.loc[wc.patient_id, 'STARcounts_RefSeq_path'],
    output:
        gene_signatures = config['results_folder'] + "/{patient_id}/05_TME/{patient_id}_gene_signatures.csv",
    threads: 1
    log: config['logs_folder'] +'/05_TME/{patient_id}_gene_signature_scores.log'
    benchmark: 'benchmarks/05_TME/{patient_id}_gene_signature_scores.log'  
    params:
        compute_scores_script = module_path + '/scripts/calculate_signature_scores.py',
        signatures_file = config['gene_signatures_file'],
        PC1_file = config['gene_signatures_PC1_file'],
        zscore_file = config['gene_signatures_zscore_file']
    shell:
        """
        (
        mkdir -p $(dirname {output.gene_signatures})

        python {params.compute_scores_script} \
                --counts_file {input} \
                --signatures_file {params.signatures_file} \
                --zscore_file {params.zscore_file} \
                --PC1_file {params.PC1_file} \
                --output_file {output.gene_signatures}
        ) > {log} 2>&1
        """