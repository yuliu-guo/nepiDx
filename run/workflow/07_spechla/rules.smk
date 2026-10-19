module_path = 'workflow/07_spechla'

rule extract_hla_read:
    input:
        bam = rules.convert_cram_to_bam_tumor_DNA.output.bam,
    output:
        hla_read_fq_1 = config["results_folder"] + "/{patient_id}/07_spechla/{patient_id}_extract_1.fq.gz",
        hla_read_fq_2 = config["results_folder"] + "/{patient_id}/07_spechla/{patient_id}_extract_2.fq.gz",
    params:
        ExtractHLAread_script = "/ref/SpecHLA/1.0.12-foss-2024a/script/ExtractHLAread.sh", # can be downloaded from SpecHLA github
        fasta_ref = config['reference']['GRCh38_fasta'],
        
    log:
        config["logs_folder"] + "/run_spechla/extract_hla_read_{patient_id}.log"
    shell:
        """
        (
        ml purge
        ml SpecHLA/1.0.12-foss-2024a
        #java -jar "$EBROOTPICARD/picard.jar"
        
        mkdir -p $(dirname {output.hla_read_fq_1}) 

        sh {params.ExtractHLAread_script} \
            -s {wildcards.patient_id} \
            -b {input.bam} \
            -r hg38 \
            -T {params.fasta_ref} \
            -o $(dirname {output.hla_read_fq_1}) 
            
        ) > {log} 2>&1

        """
        
rule run_spechla:
    input:
        hla_read_fq_1 = rules.extract_hla_read.output.hla_read_fq_1,
        hla_read_fq_2 = rules.extract_hla_read.output.hla_read_fq_2,
    output:
        freq_list = config["results_folder"] + "/{patient_id}/07_spechla/freq.list", 
    params:
        tmp = get_scratch_directory(['run_spechla']),
    log:
        config["logs_folder"] + "/run_spechla/run_spechla_{patient_id}.log"
    shell:
        """
        (
        ml purge
        ml SpecHLA/1.0.12-foss-2024a
        #java -jar "$EBROOTPICARD/picard.jar"

        spechla -n {wildcards.patient_id} \
            -1 {input.hla_read_fq_1} \
            -2 {input.hla_read_fq_2} \
            -o {params.tmp} 
        
        ls {params.tmp}/{wildcards.patient_id}/*_freq.txt > {output.freq_list}
        
        cp -a {params.tmp}/{wildcards.patient_id}/. $(dirname {output.freq_list})

        ) > {log} 2>&1

        """
        

        
rule spechla_loh:
    input:
        freq_list = rules.run_spechla.output.freq_list,
        purple_flag =rules.run_purple.output.flag,
        #hla_result = config["results_folder"] + "/{patient_id}/07_spechla/hla.result.txt", 
    output:
        loh_result = config["results_folder"] + "/{patient_id}/07_spechla/merge.hla.copy.txt", 
    params:
        tmp = get_scratch_directory(['spechla_loh']),
        spechla_loh_script = "/ref/SpecHLA/1.0.12-foss-2024a/script/cal.hla.copy.pl", # can be downloaded from SpecHLA github
        purple_results = config['results_folder'] + '/{patient_id}/07_purple/{patient_id}_tumor_DNA.purple.purity.tsv',
    log:
        config["logs_folder"] + "/run_spechla/spechla_loh_{patient_id}.log"
    shell:
        """
        (
        ml purge
        ml SpecHLA/1.0.12-foss-2024a
        #java -jar "$EBROOTPICARD/picard.jar"

        read purity ploidy < <(
            awk -F '\t' '
                NR == 1 {{
                    for (i = 1; i <= NF; i++) {{
                        if ($i == "purity") purity_col = i
                        if ($i == "ploidy") ploidy_col = i
                    }}
                }}
                NR == 2 {{
                    print $purity_col, $ploidy_col
                    exit
                }}
            ' "{params.purple_results}"
        )

        echo "Purity: $purity"
        echo "Ploidy: $ploidy"

        perl {params.spechla_loh_script} \
            -purity "$purity" \
            -ploidy "$ploidy" \
            -S {wildcards.patient_id} \
            -F {input.freq_list} \
            -T $(dirname {output.loh_result})/hla.result.txt \
            -O $(dirname {output.loh_result}) 

        ) > {log} 2>&1

        """