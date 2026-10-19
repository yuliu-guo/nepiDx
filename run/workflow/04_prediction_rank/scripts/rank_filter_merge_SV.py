import argparse
import os
import pandas as pd

def parse_arguments():
    parser = argparse.ArgumentParser(description="For SV: Filtering prediction results DataFrames from CSV files based on rank, and then merge with other summary columns")

    # CSV file arguments
    parser.add_argument(
        "-pcsv", "--prediction_csv_paths",
        nargs="+",
        required=True,
        help="Paths to the binding affinity/stability prediction CSV files, separated by space"
    )

    parser.add_argument(
        "-scsv", "--summary_csv_path",
        type=str,
        required=True,
        help="Path to the summary CSV file"
    )

    parser.add_argument(
        "-o", "--filtered_merged_output_paths",
        nargs="+",
        required=True,
        help="Paths to the output files for filtered and merged dataframes, separated by space. Need to correspond to the input prediction CSV files."
    )

    # Rank filter threshold
    parser.add_argument(
        "-r", "--rank_threshold",
        type=float,
        default=1.0,
        help="Threshold for filtering prediction CSV by rank. Default is 1.0."
    )
    return parser.parse_args()
    
def extract_tool_from_filename(path: str) -> str:
    base = os.path.basename(path)      
    stem = base.replace(".csv", "")
    tool_name = stem.split("_")[-1]    
    return tool_name   

def get_final_peptides(filtered_df_dict, original_df_dict):
    # Get all peptides that have passed the rank threshold in any of the tools
    final_peptides_tuples = set.union(*[set(zip(df['Peptide'], df['HLA_allele'], df['Identity'])) for df in filtered_df_dict.values()])

    # Get the full information for these peptides from the original dataframes
    result = {}
    for tool_name, df in original_df_dict.items():
        for peptide_tuple in final_peptides_tuples:
            peptide, hla_allele, identity = peptide_tuple
            peptide_info = df[(df['Peptide'] == peptide) & (df['HLA_allele'] == hla_allele) & (df['Identity'] == identity)]
            if tool_name not in result:
                result[tool_name] = peptide_info
            else:
                result[tool_name] = pd.concat([result[tool_name], peptide_info])

    return result

def prepare_and_filter_dataframes(prediction_csv_paths, rank_threshold, output_paths):
    # Columns to filter by in each tool
    # For now, netMHCpan rank is only taken from Rank_EL
    # TODO: We might want to include both EL and BA in the database
    tools_rank_column = {
        "netMHCpan": "Rank_EL",
        "netMHCstabpan": "%Rank_comb", #Or: "%Rank_stab", "%Rank_aff"
        "MHCflurry": "Presentation_percentile", #Or: "Affinity_percentile"
        "PRIME2": "%Rank"
    }
    output_paths_dict = {}

    filtered_df_dict = {}
    original_df_dict = {}

    # Filtering prediction dataframes
    for prediction_csv_path, output_path in zip(prediction_csv_paths, output_paths):
        df_pred = pd.read_csv(prediction_csv_path, sep = '\t')
        #tool_name = df_pred['Prediction_tool'].iloc[0]
        tool_name = extract_tool_from_filename(prediction_csv_path)
        rank_column = tools_rank_column.get(tool_name)
        if rank_column is None:
            print(f"Unknown tool name: {tool_name}. Skipping file: {prediction_csv_path}")
            continue
        df_filtered = df_pred[df_pred[rank_column] < rank_threshold]
        filtered_df_dict[tool_name] = df_filtered
        original_df_dict[tool_name] = df_pred

        # Save output paths for each tool
        output_paths_dict[tool_name] = output_path

    return filtered_df_dict, original_df_dict, output_paths_dict

def filter_and_merge_dataframes_SV(prediction_csv_paths, summary_csv_path, rank_threshold, output_paths):
    filtered_df_dict, original_df_dict, output_paths_dict = prepare_and_filter_dataframes(
        prediction_csv_paths, rank_threshold, output_paths
    )

    # Merging and saving results
    df_summary = pd.read_csv(summary_csv_path, sep = '\t')
    result_dfs = get_final_peptides(filtered_df_dict, original_df_dict) or {}
    
    # match with SV summary column names
    summary_cols = ['ID', 'sv_sorted_coord', 'sv_gene', 'svlen','svtype', 'frameshift_status','sv_gene_fusion_boolean', 'mut_kmers',
                'wt_kmers', 'novel_kmers']
    for tool_name, out_path in output_paths_dict.items():
        df_tool = result_dfs.get(tool_name)

        # minimal guard: if missing/empty, create empty df with Identity so merge + to_csv work
        if df_tool is None or df_tool.empty:
            df_tool = pd.DataFrame(columns=["Identity"])

        df_merged = df_tool.merge(
            df_summary[summary_cols],
            left_on="Identity", right_on="ID", how="left"
        ).drop(columns=["ID"])

        df_merged.to_csv(out_path, sep='\t', index=False)
        print(f"Filtered and merged data (SV) saved to {out_path}")
                
if __name__ == "__main__":
    args = parse_arguments()
    print("Processing SV data...")
    filter_and_merge_dataframes_SV(
        prediction_csv_paths=args.prediction_csv_paths,
        summary_csv_path=args.summary_csv_path,
        rank_threshold=args.rank_threshold,
        output_paths=args.filtered_merged_output_paths
    )
