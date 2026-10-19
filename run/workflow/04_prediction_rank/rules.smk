module_path = 'workflow/04_prediction_rank'



rule SNV_rank_filter_merge_netMHCpan:
    input: 
        summary_csv = rules.run_SNV_summary.output.snv_summary,
        prediction_csv = config['results_folder'] + '/{patient_id}/03_pMHC_prediction/snv_{patient_id}_netMHCpan.csv', 
    output:
        filtered_merged_output = config['results_folder'] + '/{patient_id}/04_neoepitopes/snv_{patient_id}_netMHCpan_merged_rank_filtered.csv',
    threads: 1
    log:  config['logs_folder'] +'/snv_filter_merge/{patient_id}.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} > {log} 2>&1
        ) > {log} 2>&1

        """ 

rule SNV_rank_filter_merge_netMHCstabpan:
    input: 
        summary_csv = rules.run_SNV_summary.output.snv_summary,
        prediction_csv = config['results_folder'] + '/{patient_id}/03_pMHC_prediction/snv_{patient_id}_netMHCstabpan.csv', 
    output:
        filtered_merged_output = config['results_folder'] + '/{patient_id}/04_neoepitopes/snv_{patient_id}_netMHCstabpan_merged_rank_filtered.csv',
    threads: 1
    log: config['logs_folder'] +'/snv_filter_merge/{patient_id}.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} > {log} 2>&1
        ) > {log} 2>&1

        """ 

rule SNV_rank_filter_merge_MHCflurry:
    input: 
        summary_csv = rules.run_SNV_summary.output.snv_summary,
        prediction_csv = config['results_folder'] + '/{patient_id}/03_pMHC_prediction/snv_{patient_id}_MHCflurry.csv', 
    output:
        filtered_merged_output = config['results_folder'] + '/{patient_id}/04_neoepitopes/snv_{patient_id}_MHCflurry_merged_rank_filtered.csv',
    threads: 1
    log: config['logs_folder'] +'/snv_filter_merge/{patient_id}.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} > {log} 2>&1
        ) > {log} 2>&1

        """ 

rule SNV_rank_filter_merge_PRIME2:
    input:
        summary_csv = rules.run_SNV_summary.output.snv_summary,
        prediction_csv = config['results_folder'] + '/{patient_id}/03_pMHC_prediction/snv_{patient_id}_PRIME2.csv', 
    output:
        filtered_merged_output = config['results_folder'] + '/{patient_id}/04_neoepitopes/snv_{patient_id}_PRIME2_merged_rank_filtered.csv',
    threads: 1
    log: config['logs_folder'] +'/snv_filter_merge/{patient_id}.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} > {log} 2>&1
        ) > {log} 2>&1

        """ 

rule genefusion_rank_filter_merge_netMHCpan:
    input: 
        summary_csv = rules.run_genefusion_summary.output.fusion_summary,
        prediction_csv = config['results_folder'] + '/{patient_id}/03_pMHC_prediction/fusion_{patient_id}_netMHCpan.csv',
    output:
        filtered_merged_output =  config['results_folder'] + '/{patient_id}/04_neoepitopes/fusion_{patient_id}_netMHCpan_merged_rank_filtered.csv',
    threads: 1
    log: config['logs_folder'] +'/fusion_filter_merge/{patient_id}.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} \
            -f 
        ) > {log} 2>&1

        """ 


rule genefusion_rank_filter_merge_netMHCstabpan:
    input: 
        summary_csv = rules.run_genefusion_summary.output.fusion_summary,
        prediction_csv = config['results_folder'] + '/{patient_id}/03_pMHC_prediction/fusion_{patient_id}_netMHCstabpan.csv',
    output:
        filtered_merged_output =  config['results_folder'] + '/{patient_id}/04_neoepitopes/fusion_{patient_id}_netMHCstabpan_merged_rank_filtered.csv',
    threads: 1
    log: config['logs_folder'] +'/fusion_filter_merge/{patient_id}.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} \
            -f 
        ) > {log} 2>&1

        """ 


rule genefusion_rank_filter_merge_MHCflurry:
    input: 
        summary_csv = rules.run_genefusion_summary.output.fusion_summary,
        prediction_csv = config['results_folder'] + '/{patient_id}/03_pMHC_prediction/fusion_{patient_id}_MHCflurry.csv',
    output:
        filtered_merged_output =  config['results_folder'] + '/{patient_id}/04_neoepitopes/fusion_{patient_id}_MHCflurry_merged_rank_filtered.csv',
    threads: 1
    log: config['logs_folder'] +'/fusion_filter_merge/{patient_id}.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} \
            -f 
        ) > {log} 2>&1

        """ 



rule genefusion_rank_filter_merge_PRIME2:
    input: 
        summary_csv = rules.run_genefusion_summary.output.fusion_summary,
        prediction_csv = config['results_folder'] + '/{patient_id}/03_pMHC_prediction/fusion_{patient_id}_PRIME2.csv',
    output:
        filtered_merged_output =  config['results_folder'] + '/{patient_id}/04_neoepitopes/fusion_{patient_id}_PRIME2_merged_rank_filtered.csv',
    threads: 1
    log: config['logs_folder'] +'/fusion_filter_merge/{patient_id}.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} \
            -f 
        ) > {log} 2>&1

        """ 



rule sv_rank_filter_merge_netMHCpan:
    input: 
        #unpack(sv_rank_filter_merge_netMHCpan)
        summary_csv = rules.run_sv_summary.output.out,
        prediction_csv = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/sv_{patient_id}_netMHCpan.csv'),
    output:
        filtered_merged_output =  (config['results_folder'] + '/{patient_id}/04_neoepitopes/sv_{patient_id}_netMHCpan_merged_rank_filtered.csv'),
    threads: 1
    log: config['logs_folder'] +'/sv_filter_merge/{patient_id}_netMHCpan.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge_SV.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} \
        ) > {log} 2>&1

        """ 

rule sv_rank_filter_merge_netMHCstabpan:
    input: 
        summary_csv = rules.run_sv_summary.output.out,
        prediction_csv = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/sv_{patient_id}_netMHCstabpan.csv'),
    output:
        filtered_merged_output =  (config['results_folder'] + '/{patient_id}/04_neoepitopes/sv_{patient_id}_netMHCstabpan_merged_rank_filtered.csv'),
    threads: 1
    log: config['logs_folder'] +'/sv_filter_merge/{patient_id}_netMHCstabpan.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge_SV.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} \
        ) > {log} 2>&1

        """ 

rule sv_rank_filter_merge_MHCflurry:
    input: 
        summary_csv = rules.run_sv_summary.output.out,
        prediction_csv = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/sv_{patient_id}_MHCflurry.csv'),
    output:
        filtered_merged_output =  (config['results_folder'] + '/{patient_id}/04_neoepitopes/sv_{patient_id}_MHCflurry_merged_rank_filtered.csv'),
    threads: 1
    log: config['logs_folder'] +'/sv_filter_merge/{patient_id}_MHCflurry.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge_SV.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} \
        ) > {log} 2>&1

        """ 

rule sv_rank_filter_merge_PRIME2:
    input: 
        summary_csv = rules.run_sv_summary.output.out,
        prediction_csv = (config['results_folder'] + '/{patient_id}/03_pMHC_prediction/sv_{patient_id}_PRIME2.csv'),
    output:
        filtered_merged_output =  (config['results_folder'] + '/{patient_id}/04_neoepitopes/sv_{patient_id}_PRIME2_merged_rank_filtered.csv'),
    threads: 1
    log: config['logs_folder'] +'/sv_filter_merge/{patient_id}_PRIME2.log'
    params:
        rank_threshold = 1,
        scripts =  module_path + '/scripts/rank_filter_merge_SV.py',
    shell:
        """
        (
        python {params.scripts} \
            -scsv {input.summary_csv} \
            -pcsv {input.prediction_csv} \
            -o {output.filtered_merged_output} \
            -r {params.rank_threshold} \
        ) > {log} 2>&1

        """ 