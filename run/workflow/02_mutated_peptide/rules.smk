module_path = 'workflow/02_mutated_peptide'

rule peptide_similarity_SNV:
    input:
        snv_summary = rules.run_SNV_summary.output.snv_summary,
    output:
        snv_peptides_sim =  (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_snv_similarity.csv'),
    log: config['logs_folder'] +'/calculate_self_similarity_snv/{patient_id}_similarity.log'
    benchmark: 'benchmarks/calculate_self_similarity_snv/{patient_id}_similarity.log'  
    params:
        scripts = module_path + '/scripts/expand_and_match.py',
        seq_col = 'flanked_mutated_protein_sequence',
        keep_col = 'ID SYMBOL RNA_ALT_N',
        #k = config['params']['k_mer'],
        k = 9,
        db_path = config['human_k9'],
        blosum62_2 = module_path + '/scripts/blosum62-2.qij',
        k1_cache = module_path + '/scripts/blosum62-2.qij.beta_0.125.k1_cache.npz',
    shell:
        """
        (
	python {params.scripts} \
            --input-tsv {input.snv_summary} \
            --output {output.snv_peptides_sim} \
            --seq-col {params.seq_col} \
            --keep-cols {params.keep_col} \
            --kmer {params.k} \
            --db {params.db_path} \
            --blos {params.blosum62_2} \
            --cache-path {params.k1_cache}
	) > {log} 2>&1
        """

rule collect_mutated_peptide_SNV:
    input:
        snv_peptides_sim = rules.peptide_similarity_SNV.output.snv_peptides_sim,
    output:
        #snv_flanked_fasta =  (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_snv_flanked.fasta'),
        snv_nmer_fasta =  (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_snv_nmer.fasta'),
    threads: 1
    log: config['logs_folder'] +'/filter_self_peptide_snv/{patient_id}.log'
    benchmark: 'benchmarks/filter_self_peptide_snv/{patient_id}.log'  
    params:
        scripts =  module_path + '/scripts/filter_snv_self_peptides.py',
        min_alt_base = config['params']['min_alt_base_snv'],
        max_hamming_similarity = config['params']['max_hamming_similarity'],
        max_blosum62_2_similarity = config['params']['max_blosum62_2_similarity'],
    shell:
        """
        (
        python {params.scripts} \
            --input_path {input.snv_peptides_sim} \
            --nmer_output_path {output.snv_nmer_fasta} \
            --min_alt_read_base {params.min_alt_base} \
            --max_hamming_similarity {params.max_hamming_similarity} \
            --max_blosum62_2_similarity {params.max_blosum62_2_similarity}
        ) > {log} 2>&1

        """

rule peptide_similarity_fusion:
    input:
        fusion_summary = rules.run_genefusion_summary.output.fusion_summary,
    output:
        fusion_peptides_sim =  (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_fusion_similarity.csv'),
    log: config['logs_folder'] +'/calculate_self_similarity_fusion/{patient_id}_similarity.log'
    benchmark: 'benchmarks/calculate_self_similarity_fusion/{patient_id}_similarity.log'  
    params:
        scripts = module_path + '/scripts/expand_and_match.py',
        seq_col = 'fusion_amino_acid_sequence',
        keep_col = 'ID gene_name fusion_read_support fusion_expression',
        #k = config['params']['k_mer'],
        k = 9,
        db_path = config['human_k9'],
        #db_h_k9 = module_path + '/scripts/human_k9.sqlite',
        blosum62_2 = module_path + '/scripts/blosum62-2.qij',
        k1_cache = module_path + '/scripts/blosum62-2.qij.beta_0.125.k1_cache.npz',
    shell:
        """
        (python {params.scripts} \
            --input-tsv {input.fusion_summary} \
            --output {output.fusion_peptides_sim} \
            --seq-col {params.seq_col} \
            --keep-cols {params.keep_col} \
            --kmer {params.k} \
            --db {params.db_path} \
            --blos {params.blosum62_2} \
            --cache-path {params.k1_cache}
	) > {log} 2>&1
        """


rule collect_mutated_peptide_fusion:
    input:
        fusion_peptides_sim = rules.peptide_similarity_fusion.output.fusion_peptides_sim,
    output:
        #fusion_flanked_fasta =  (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_fusion_flanked.fasta'),
        fusion_nmer_fasta =  (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_fusion_nmer.fasta'),
    threads: 1
    log: config['logs_folder'] +'/filter_self_peptide_fusion/{patient_id}.log'
    benchmark: 'benchmarks/filter_self_peptide_fusion/{patient_id}.log'  
    params:
        scripts =  module_path + '/scripts/filter_fusion_self_peptides.py',
        min_fusion_read_support = config['params']['min_fusion_read_support'],
        max_hamming_similarity = config['params']['max_hamming_similarity'],
        max_blosum62_2_similarity = config['params']['max_blosum62_2_similarity'],
    shell:
        """
        (
        python {params.scripts} \
            --input_path {input.fusion_peptides_sim} \
            --nmer_output_path {output.fusion_nmer_fasta} \
            --min_fusion_read_support {params.min_fusion_read_support} \
            --max_hamming_similarity {params.max_hamming_similarity} \
            --max_blosum62_2_similarity {params.max_blosum62_2_similarity}
        ) > {log} 2>&1
        """

rule peptide_similarity_sv:
    input:
        sv_summary = rules.run_sv_summary.output.out,
    output:
        sv_peptides_sim =  (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_sv_similarity.csv'),
    log: config['logs_folder'] +'/calculate_self_similarity_sv/{patient_id}_similarity.log'
    benchmark: 'benchmarks/calculate_self_similarity_sv/{patient_id}_similarity.log'  
    params:
        scripts = module_path + '/scripts/expand_and_match_sv.py',
        seq_col = 'novel_kmers_sequence',
        keep_col = 'ID sv_gene',
        #k = config['params']['k_mer'],
        k = 9,
        db_path = config['human_k9'],
        #db_h_k9 = module_path + '/scripts/human_k9.sqlite',
        blosum62_2 = module_path + '/scripts/blosum62-2.qij',
        k1_cache = module_path + '/scripts/blosum62-2.qij.beta_0.125.k1_cache.npz',
    shell:
        """
        (
	python {params.scripts} \
            --input-tsv {input.sv_summary} \
            --output {output.sv_peptides_sim} \
            --seq-col {params.seq_col} \
            --keep-cols {params.keep_col} \
            --kmer {params.k} \
            --db {params.db_path} \
            --blos {params.blosum62_2} \
            --cache-path {params.k1_cache}
	) > {log} 2>&1
        """

        
rule collect_mutated_peptide_sv:
    input:
        #sv_summary = rules.run_sv_summary.output.out,
        sv_peptides_sim = rules.peptide_similarity_sv.output.sv_peptides_sim,
        #sv_summary = (config['results_folder'] + '/{patient_id}/00_variants_summary/{patient_id}_sv_summary.csv'),
    output:
        sv_nmer_fasta = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_sv_nmer.fasta'),
    threads: 1
    log: config['logs_folder'] +'/filter_self_peptide_sv/{patient_id}.log'
    benchmark: 'benchmarks/filter_self_peptide_sv/{patient_id}.log'  
    params:
        scripts =  module_path + '/scripts/filter_sv_self_peptides.py',
        max_hamming_similarity = config['params']['max_hamming_similarity'],
        max_blosum62_2_similarity = config['params']['max_blosum62_2_similarity'],
    shell:
        """
        (
        python {params.scripts} \
            --input_path {input.sv_peptides_sim} \
            --nmer_output_path {output.sv_nmer_fasta} \
            --max_hamming_similarity {params.max_hamming_similarity} \
            --max_blosum62_2_similarity {params.max_blosum62_2_similarity}
        ) > {log} 2>&1
        """
