module_path = 'workflow/09_diann'
import csv
from pathlib import Path


def read_tsv_mapping(path, key_column, value_column):
    """Read a two-column mapping from a tab-separated file."""
    mapping = {}

    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required_columns = {key_column, value_column}
        missing_columns = required_columns - set(reader.fieldnames or [])

        if missing_columns:
            raise ValueError(
                f"Missing columns in {path}: {sorted(missing_columns)}"
            )

        for row_number, row in enumerate(reader, start=2):
            key = row[key_column].strip()
            value = row[value_column].strip()

            # Ignore rows with a missing patient, plate, or path.
            if not key or not value:
                continue

            if key in mapping and mapping[key] != value:
                raise ValueError(
                    f"Conflicting values for {key!r} in {path},"
                    f"around line {row_number}: "
                    f"{mapping[key]!r} versus {value!r}"
                )

            mapping[key] = value

    return mapping


PATIENT_TO_PLATE = read_tsv_mapping(
    "/ref/MUT_fasta_proteomic/protarget_65_patient_plate.tsv",
    key_column="patient_id",
    value_column="plate_nr",
)

PLATE_TO_MZML = read_tsv_mapping(
    "/ref/MUT_fasta_proteomic/protarget_65_plate_mzml.tsv",
    key_column="plate_nr",
    value_column="mzml_path",
)


def build_patient_mzml_mapping():
    """
    Include only patients that have:
      a plate number,
      and an mzML path for that plate.
    """
    patient_to_mzml = {}

    for patient_id, plate_nr in PATIENT_TO_PLATE.items():
        mzml_path = PLATE_TO_MZML.get(plate_nr)

        if mzml_path:
            patient_to_mzml[patient_id] = mzml_path

    return patient_to_mzml


PATIENT_TO_MZML = build_patient_mzml_mapping()

# Only these patients should be used when expanding DIA-NN targets.
DIANN_PATIENTS = sorted(PATIENT_TO_MZML)


def get_mzml_from_patient(wildcards):
    """Snakemake input function for the DIA-NN mzML input."""
    patient_id = wildcards.patient_id

    try:
        return PATIENT_TO_MZML[patient_id]
    except KeyError as error:
        raise ValueError(
            f"No mzML file is available for patient_id={patient_id!r}. "
            "The patient may have no plate number, or the plate may have "
            "no matching mzML file."
        ) from error

rule run_diann:
    input:
        mut_fasta = (config['results_folder'] + '/fasta_MUT/UP000005640_9606_HERV_orfs_aa_sequences_v5_header_modified_{patient_id}_MUT.fasta'),
        mzml = get_mzml_from_patient,
    output:
        report= (config['results_folder'] + '/{patient_id}/09_diann/{patient_id}.stats.tsv'),
        #outlib=(config['results_folder'] + '/{patient_id}/09_diann/{patient_id}_outlib.tsv'),
        #outdir=(config['results_folder'] + '/{patient_id}/09_diann/'),
    threads:20
    log: config['logs_folder'] +'/run_diann/{patient_id}.log'
    benchmark: 'benchmarks/run_diann/{patient_id}.log'
    params:
        tmp = get_scratch_directory(['diann']),
    shell:
        """
        (
        ml DIA-NN/2.3.2-gfbf-2024a

        mkdir -p $(dirname {output.report}) &&
        mkdir -p {params.tmp} &&

        diann-2.3.2 \
            --f {input.mzml} \
            --out {params.tmp}/{wildcards.patient_id}.tsv \
            --out-lib {params.tmp}/{wildcards.patient_id}_outlib.tsv \
            --fasta {input.mut_fasta} \
            --temp {params.tmp} \
            --gen-spec-lib \
            --predictor \
            --verbose 1 \
            --qvalue 0.01 \
            --matrices \
            --xic \
            --threads 20 \
            --fasta-search \
            --min-pep-len 7 \
            --max-pep-len 52 \
            --min-pr-mz 300 \
            --max-pr-mz 1800 \
            --min-pr-charge 2 \
            --max-pr-charge 7 \
            --cut 'K*,R*' \
            --missed-cleavages 2 \
            --var-mods 2 \
            --var-mod 'UniMod:35,15.994915,M' \
            --var-mod 'UniMod:1,42.010565,*n' \
            --mass-acc 15 \
            --mass-acc-ms1 15 \
            --use-quant \
            --peptidoforms \
            --reanalyse \
            --relaxed-prot-inf \
            --smart-profiling
            
            
        cp -a {params.tmp}/. $(dirname {output.report})

        ) > {log} 2>&1
        """



rule find_mut_peptide:
    input:
        diann_report = rules.run_diann.output.report,
        #diann_report=(config['results_folder'] + '/{patient_id}/09_diann/{patient_id}.parquet'),
    output:
        flag = (config['results_folder'] + '/{patient_id}/09_diann/mutation_specific_results/all_tagged_peptides_classification.tsv'),
    threads:20
    log: config['logs_folder'] +'/run_diann/{patient_id}_find_mut_peptide.log'
    benchmark: 'benchmarks/run_diann/{patient_id}_find_mut_peptide.log'
    params:
        script = module_path + '/scripts/find_mutation_specific_peptides_from_diann.py',
    shell:
        """
        (
        ml Arrow/22.0.0-gfbf-2025a
        
        python {params.script} $(dirname {input.diann_report})

        ) > {log} 2>&1
        """
        
        
rule find_target_peptide:
    input:
        diann_report = rules.run_diann.output.report,
    output:
        flag = (config['results_folder'] + '/{patient_id}/09_diann/target_specific_results/all_tagged_peptides_classification.tsv'),
    threads:20
    log: config['logs_folder'] +'/run_diann/{patient_id}_find_target_peptide.log'
    benchmark: 'benchmarks/run_diann/{patient_id}_find_target_peptide.log'
    params:
        script = module_path + '/scripts/find_target_specific_peptides_from_diann.py',
    shell:
        """
        (
        ml Arrow/22.0.0-gfbf-2025a
        
        python {params.script} $(dirname {input.diann_report})

        ) > {log} 2>&1
        """
        
#rule find_target_mut_peptide_all:
#    input:
#        mut_flag = 
#        target_flag = 
#    output: 
        
rule gather_and_index_ms_peptide:
    input:
        mut_peptide = (config['results_folder'] + '/{patient_id}/09_diann/mutation_specific_results/mutation_specific_peptides.tsv'),
        target_peptide = (config['results_folder'] + '/{patient_id}/09_diann/target_specific_results/target_specific_peptides.tsv'),
    output:
        summary_file = (config['results_folder'] + '/{patient_id}/09_diann/{patient_id}_summary_specific_peptides.tsv'),
    log: config['logs_folder'] +'/run_diann/{patient_id}_gather_and_index_ms_peptide.log'
    benchmark: 'benchmarks/run_diann/{patient_id}_gather_and_index_ms_peptide.log'
    params:
        script = module_path + '/scripts/gather_and_index_ms_peptide.py',
    shell:
        """
        (

        python {params.script} \
            --target_specific_file {input.target_peptide} \
            --mut_specific_file {input.mut_peptide} \
            --output_file {output.summary_file}

        ) > {log} 2>&1
        """
        
rule find_all_9mer_from_ms_peptide:
    input:
        peptides = (config['results_folder'] + '/{patient_id}/09_diann/{patient_id}_summary_specific_peptides.tsv'),
    output:
        fasta_9mer = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_ms9mer_nmer.fasta'),
    log: config['logs_folder'] +'/run_diann/{patient_id}_find_all_9mer_from_ms_peptide.log'
    benchmark: 'benchmarks/run_diann/{patient_id}_find_all_9mer_from_ms_peptide.log'
    params:
        script = module_path + '/scripts/generate_nmers.py',
        nmer_length = 9
    shell:
        """
        (

        python {params.script} \
            --input_file {input.peptides} \
            --nmer_length {params.nmer_length} \
            --output_file {output.fasta_9mer}

        ) > {log} 2>&1
        """
        
rule find_all_8mer_from_ms_peptide:
    input:
        peptides = (config['results_folder'] + '/{patient_id}/09_diann/{patient_id}_summary_specific_peptides.tsv'),
    output:
        fasta_8mer = (config['results_folder'] + '/{patient_id}/02_mutated_peptide/{patient_id}_ms8mer_nmer.fasta'),
    log: config['logs_folder'] +'/run_diann/{patient_id}_find_all_8mer_from_ms_peptide.log'
    benchmark: 'benchmarks/run_diann/{patient_id}_find_all_8mer_from_ms_peptide.log'
    params:
        script = module_path + '/scripts/generate_nmers.py',
        nmer_length = 8
    shell:
        """
        (

        python {params.script} \
            --input_file {input.peptides} \
            --nmer_length {params.nmer_length} \
            --output_file {output.fasta_8mer}

        ) > {log} 2>&1
        """