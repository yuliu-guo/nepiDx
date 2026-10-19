module_path = 'workflow/initial_setup.smk'

import re
import os
import pathlib
import pandas as pd
from typing import List, Iterable, Optional

def get_scratch_directory(subdirs: List[str], force_workdir_scratch: bool=False) -> str:
    '''
    Auxiliary function that returns the scratch directory location
    for a given rule.

    Arguments
    ---------
    subdirs: List[str]
        List of subdirectory names.
    force_workdir_scratch: bool
        If true, force the scratch directory inside of the
        working directory.
    '''
    jobid = os.environ.get('PBS_JOBID', None)
    basedir = None
    if not jobid or force_workdir_scratch:
        basedir = Path(pipeline_workdir).joinpath('scratch')
    else:
        basedir = Path('/scratch').joinpath(jobid)
    if subdirs:
        for element in subdirs:
            basedir = basedir.joinpath(element)
    return basedir.as_posix()


def parse_file_path(file_path):
    filepath = Path(file_path)
    filename = filepath.stem  # Get filename without extension
    
    # Extract sample_id from the filename
    sample_id = re.sub(r'_L[0-9]{3}', '', filename)
    sample_id = re.sub(r'_[cs]nv_.*$', '', sample_id)
    sample_id = re.sub(r'_mqc$', '', sample_id)

    return sample_id  

## Flatten a list of lists
def flatten(list_of_lists: List) -> List:
	result = list()
	for i in list_of_lists:
		if isinstance(i, list):
			result.extend(flatten(i))
		else:
			result.append(str(i))
	return result

## Get bind string for singularity exec
def get_bind_string(input, output, extra, scratch_dir, singularity_internal_dir="", add_singularity_dir=False):
  to_process = []
  to_process.extend(flatten(input))
  to_process.extend(flatten(output))
  to_process.extend(flatten(extra))
     
  dirs = set()
  if scratch_dir and str(scratch_dir).strip():
    dirs.add(scratch_dir)
    
  for element in to_process:
    if not element:
      continue

    path = Path(element).absolute()
    
    # If input is a symlink, also bind where it points
    if path.is_symlink():
      try:
        target = path.resolve(strict=True)
        if os.path.isdir(target):
          dirs.add(str(target))
        else:
          dirs.add(target.parent.as_posix())
      except Exception:
        pass
            
    if os.path.isdir(path):
      dirs.add(str(path))
    else:
      dirs.add(path.parent.as_posix())
    
  result = ' '.join([f'--bind {path}' for path in sorted(dirs)])
    
  if add_singularity_dir:
    if scratch_dir and str(scratch_dir).strip() and singularity_internal_dir:
      result = result + ' --bind ' + scratch_dir + ':' + singularity_internal_dir
    
  return result
  
def get_inputs_with_quotation(input_list):
  input_paths = flatten(input_list)
  result = []
  for path in input_paths:
    result.append(f'"{path}"')
    
  return ' '.join(result)

def prune_rule_if_empty_input(*rules):
    
    for r in rules:
        if not isinstance(r, snakemake.rules.RuleProxy):
            raise ValueError('The received object is not a snakemake rule.')
        if not r.input:
            r._output = []

def hla_allele_path(wc):
    # Raw value from the dataframe
    val = df.loc[wc.patient_id, "HLA_allele"]

    # Case 1: real NaN (missing in table)
    if pd.isna(val):
        return []

    # Normalize to clean string
    s = str(val).strip()

    # Case 2: empty string or sentinel values
    if s == "" or s.lower() in ("none", "nan"):
        return []

    # Case 3: real path
    return s

def check_file_exist(wc, column:str):
    #check if the file exist for certain filetype (HLA_allele, RNA_cram_path, ects)
    val = df.loc[wc.patient_id, column]

    if pd.isna(val):
        return []
    s = str(val).strip()

    if s =="" or s.lower() in ("none", "nan"):
        return []
    
    return s


def get_patient_id_for_rule_all(
    df: pd.DataFrame,
    option1: Optional[Iterable[str]] = None,
    option2: str = "summary",
) -> pd.Index:
    """
    Return df.index entries (patient IDs) that satisfy required-column rules.

    option1: iterable of rule types among {"sv","snv","fusion","expression","mhchammer"} or None -> all
    option2: "summary" (default) or "prediction"
    """

    # Normalize inputs
    if option2 is None:
        option2 = "summary"
    option2 = str(option2).strip().lower()
    if option2 not in {"summary", "prediction"}:
        raise ValueError("option2 must be 'summary' or 'prediction'")

    all_rules = ["sv", "snv", "fusion", "hla_tumor_rna", "hla_tumor_dna", "hla_normal_dna", "expression", "isovar", "mhchammer", "lohhla", "spechla", "database"]
    if option1 is None:
        rules = all_rules
    else:
        rules = [str(r).strip().lower() for r in option1]
        unknown = sorted(set(rules) - set(all_rules))
        if unknown:
            raise ValueError(f"Unknown rule(s) in option1: {unknown}. Allowed: {all_rules}")

    # Required columns by (rule, mode)
    required_cols = {
        ("snv", "prediction"): ["snv_somatic_annot_vcf_path", "Paired_tumor_DNA_cram_path"],
        ("fusion", "prediction"): ["gene_fusion_arriba", "Paired_tumor_DNA_cram_path"],
        ("sv", "prediction"): ["sv_somatic_raw_manta", "Paired_tumor_DNA_cram_path"],

        ("snv", "summary"): ["snv_somatic_annot_vcf_path"],
        ("fusion", "summary"): ["gene_fusion_arriba"],
        ("sv", "summary"): ["sv_somatic_raw_manta"],
        
        # check if RNA and snv exist for isovar
        ("isovar", "summary"): ["snv_somatic_annot_vcf_path" ,"RNA_cram_path"],

        # HLA-typing
        ("hla_tumor_rna", "summary"): ["Paired_tumor_DNA_cram_path"],
        ("hla_tumor_dna", "summary"): ["Paired_tumor_DNA_cram_path"],
        ("hla_normal_dna", "summary"): ["Paired_normal_DNA_cram_path"],

        # expression: summary (default)
        ("expression", "summary"): ["STARcounts_RefSeq_path"],

        # mhchammer: summary
        ("mhchammer", "summary"): [
            "Paired_normal_DNA_cram_path",
            "RNA_cram_path",
            "Paired_tumor_DNA_cram_path",
        ],
        # lohhla: summary
        ("lohhla", "summary"): [
            "Paired_normal_DNA_cram_path",
            "RNA_cram_path",
            "Paired_tumor_DNA_cram_path",
        ],
        
        # spechla: summary
        ("spechla", "summary"): [
            "Paired_normal_DNA_cram_path",
            "Paired_tumor_DNA_cram_path",
        ],
        
        # database: summary
        ("database", "summary"): [
            "Paired_normal_DNA_cram_path",
            "RNA_cram_path",
            "Paired_tumor_DNA_cram_path",
            "snv_somatic_annot_vcf_path",
            "gene_fusion_arriba",
            "sv_somatic_raw_manta",
            "STARcounts_RefSeq_path"
        ],
    }

    # Helper: treat NaN / None / "" / whitespace as missing
    def _present_mask(series: pd.Series) -> pd.Series:
        s = series.astype("string").str.strip()
        return s.notna() & ~s.str.lower().isin({"", "na", "none", "nan"})

    matched_ids: Set = set()

    for rule in rules:
        key = (rule, option2)

        # Special default behavior for expression when option2 isn't provided: already default "summary"
        if key not in required_cols:
            # If user asked for expression or mhchammer with prediction, spec doesn't define it -> skip
            # (or raise if you prefer strict behavior)
            continue

        cols = required_cols[key]

        missing_cols = [c for c in cols if c not in df.columns]
        if missing_cols:
            raise KeyError(f"Missing required column(s) for {rule}/{option2}: {missing_cols}")

        mask = pd.Series(True, index=df.index)
        for c in cols:
            mask &= _present_mask(df[c])

        matched_ids.update(df.index[mask].tolist())

    # Return as an Index (unique, sorted for reproducibility)
    try:
        return pd.Index(sorted(matched_ids))
    except TypeError:
        # If index values aren't sortable (mixed types), return unsorted but unique
        return pd.Index(list(matched_ids))


