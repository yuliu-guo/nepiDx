import pandas as pd
import numpy as np
import argparse
import os

def housekeeping_normalization(counts_df, housekeeping_genes):
    """
    Normalize gene expression counts using housekeeping genes.
    
    Parameters:
    counts_df (pd.DataFrame): DataFrame containing gene expression counts with 'gene' as the first column.
    housekeeping_genes (list): List of housekeeping genes to use for normalization.
    
    Returns:
    pd.DataFrame: Normalized gene expression counts.
    """
    normalized_counts = counts_df.copy()
    
    # Log-transform all counts (adding a small constant to avoid log(0))
    normalized_counts['log_counts'] = np.log2(normalized_counts['STDIN'] + 1)
    
    # Filter for housekeeping genes
    hk_counts = normalized_counts[normalized_counts['gene'].isin(housekeeping_genes)]
    
    if hk_counts.empty:
        raise ValueError("No housekeeping genes found in counts_df.")
    
    # Calculate the mean expression of housekeeping genes
    hk_mean = hk_counts['log_counts'].mean()
    
    # Normalize the counts by dividing by the mean of housekeeping genes
    normalized_counts['normalized_counts'] = normalized_counts['log_counts'] - hk_mean
    
    return normalized_counts[['gene', 'normalized_counts']]
    
def zscore_transformation(counts_df, zscore_df):
    """
    Compute z-score normalized expression values for a sample using background statistics.
    
    Parameters:
    counts_df (pd.DataFrame): DataFrame containing normalized expression values for one sample.
    zscore_df (pd.DataFrame): DataFrame containing background statistics per gene.

    Returns:
    pd.DataFrame: DataFrame with zscore transformed normalized counts.
    """
    # Merge the two dataframes on gene ID
    merged = counts_df.merge(
        zscore_df,
        left_on='gene',
        right_on='gene_id',
        how='left'
    )
    
    # Compute z-score
    merged['zscore'] = (
        (merged['normalized_counts'] - merged['gene_mean']) /
        merged['gene_std']
    )
    
    # Handle division by zero or missing std
    merged.loc[(merged['gene_std'] == 0) | (merged['gene_std'].isna()), 'zscore'] = np.nan
    
    # Return only relevant columns
    return merged[['gene', 'zscore']]
    
def PC1_calculation(zscore_counts_df, PC1_df, gene_list, signature_name):
    """
    Calculate the PC1 score for a given gene signature using z-score transformed data.

    This function computes the first principal component (PC1) score for a sample by
    taking a weighted sum of z-score normalized gene expression values. The weights
    (PC1 coefficients) are provided from background samples.

    Parameters:
    zscore_counts_df (pd.DataFrame): DataFrame containing z-score transformed expression values.
    PC1_df (pd.DataFrame): DataFrame containing PC1 coefficients per gene for different signatures.
    gene_list (list): List of genes included in the gene signature.
    signature_name (str): Name of the gene signature. Must match a column in PC1_df.

    Returns:
    float: The PC1 score for the given sample and gene signature.
    """

    # Filter z-score data to genes in the gene set
    z_df = zscore_counts_df[zscore_counts_df['gene'].isin(gene_list)]

    # Filter PC1 coefficients to genes in the gene set
    pc1_sub = PC1_df[PC1_df['Geneid'].isin(gene_list)][['Geneid', signature_name]]

    # Merge on gene IDs
    merged = z_df.merge(
        pc1_sub,
        left_on='gene',
        right_on='Geneid',
        how='inner'
    )

    # Compute weighted sum (dot product)
    merged['weighted'] = merged['zscore'] * merged[signature_name]

    pc1_score = merged['weighted'].sum()

    return pc1_score

def T_cell_inflamed_GEP_calculation(counts_df):
    """
    Calculate the T cell-inflamed Gene Expression Profile (GEP) score
    as a weighted sum of normalized gene expression values.

    Parameters:
    counts_df (pd.DataFrame): DataFrame containing normalized expression values.

    Returns:
    float: Weighted sum GEP score.
    """
    
    weights_dict = {
      "gene": ["CCL5", "CD27", "CD274", "CD276","CD8A","CMKLR1","CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"],
      "weight": [0.008346, 0.072293, 0.042853, -0.0239, 0.031021, 0.151253, 0.074135, 0.004313, 0.020091, 0.058806, 0.07175, 0.060679, 0.123895, 0.075524, 0.003734, 0.032999, 0.250229, 0.084767]
    }
    weights_df = pd.DataFrame(weights_dict)
    
    tcell_df = counts_df.merge(
        weights_df,
        left_on='gene',
        right_on='gene',
        how='inner'
    )
    
    gep_score = (tcell_df['normalized_counts'] * tcell_df['weight']).sum()
    
    return gep_score
    
def geometric_mean_log2_normalized(counts_df, gene_list):
    """
    Compute the geometric mean of housekeeping-normalized expression values
    for a given gene set.
    Assumes 'normalized_counts' are in log2 space.

    Parameters:
    counts_df (pd.DataFrame): DataFrame containing normalized expression values.
    gene_list (list): List of genes in the signature

    Returns:
    float: Geometric mean in linear space (or NaN if no valid genes)
    """
    subset = counts_df[counts_df['gene'].isin(gene_list)]['normalized_counts']
    subset = subset.dropna()

    if len(subset) == 0:
        return np.nan

    return 2 ** subset.mean()

def calculate_signature_scores(counts_df, signatures_df, zscore_df, PC1_df):
    """
    Calculate gene signature scores based on normalized gene expression counts.
    
    Parameters:
    counts_df (pd.DataFrame): DataFrame containing normalized gene expression counts with 'gene' and 'normalized_counts' columns.
    signatures_df (pd.DataFrame): DataFrame containing gene signatures with 'Gene' and 'signatureID' columns.
    
    Returns:
    pd.DataFrame: DataFrame containing gene signature scores.
    """
    signature_scores = {'score': [], 'method': [], 'value': []}
    zscore_transformed = zscore_transformation(counts_df, zscore_df)
    
    for signature in signatures_df['signatureID'].unique():
        genes_in_signature = signatures_df[signatures_df['signatureID'] == signature]['Gene'].tolist()
        #print("Columns:", signatures_df.columns.tolist())
        #print(signatures_df.head())
        if 'Method' not in signatures_df.columns:
            raise ValueError(f"'Method' column not found. Available columns: {signatures_df.columns.tolist()}")
        
        method = signatures_df[signatures_df['signatureID'] == signature]['Method'].iloc[0]

        signature_scores['score'].append(signature)
        signature_scores['method'].append(method)
        
        if method == "Simple average":
            signature_scores['value'].append(counts_df[counts_df['gene'].isin(genes_in_signature)]['normalized_counts'].mean())
        elif method == "Average zscore":
            signature_scores['value'].append(zscore_transformed[zscore_transformed['gene'].isin(genes_in_signature)]['zscore'].mean())
        elif method == "PC1":
            signature_scores['value'].append(PC1_calculation(zscore_transformed, PC1_df, genes_in_signature, signature))
        elif method == "T cell GEP":
            signature_scores['value'].append(T_cell_inflamed_GEP_calculation(counts_df))
        elif method == "Geometric mean":
            signature_scores['value'].append(geometric_mean_log2_normalized(counts_df, genes_in_signature))
        else:
            raise ValueError(f"Unknown method: {method}")
            
    output_df = pd.DataFrame(signature_scores)
    
    return output_df

def main():
    parser = argparse.ArgumentParser(description='Calculate gene signature scores from gene expression counts.')
    parser.add_argument('--counts_file', type=str, required=True, help='Path to the gene expression counts file (tab-delimited).')
    parser.add_argument('--signatures_file', type=str, required=True, help='Path to the gene signatures file with gene lists and methods (tab-delimited).')
    parser.add_argument('--zscore_file', type=str, required=True, help='Path to the mean and standard deviations of expression per gene (tab-delimited).')
    parser.add_argument('--PC1_file', type=str, required=True, help='Path to the PC1 coefficients for gene sets with PC1 method (tab-delimited).')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output file for gene signature scores (tab-delimited).')
    
    args = parser.parse_args()
    
    # Load the data
    counts_df = pd.read_csv(args.counts_file, sep='\t')
    signatures_df = pd.read_csv(args.signatures_file, sep='\t')
    zscore_df = pd.read_csv(args.zscore_file, sep='\t')
    PC1_df = pd.read_csv(args.PC1_file, sep='\t')

    # Normalize accepted counts formats:
    #   1) gene \t STDIN
    #   2) gene \t unstrand \t strand1 \t strand2
    if counts_df.shape[1] == 2:
        counts_df.columns = ['gene', 'STDIN']
    elif counts_df.shape[1] == 4:
        counts_df.columns = ['gene', 'STDIN', 'strand1', 'strand2']
    else:
        raise ValueError(
            f"counts_file must have either 2 or 4 columns, found {counts_df.shape[1]} columns."
        )
    
    mask = counts_df['gene'].astype(str).str.startswith('N_')
    counts_df = counts_df[mask == False].copy()
    
    # Normalize the gene expression counts using housekeeping genes
    housekeeping_genes = ["STK11IP", "ZBTB34", "OAZ1", "POLR2A", "G6PD","ABCF1", "UBB", "TBP", "SDHA"]
    normalized_counts_df = housekeeping_normalization(counts_df, housekeeping_genes)
    
    # Calculate gene signature scores
    output_df = calculate_signature_scores(normalized_counts_df, signatures_df, zscore_df, PC1_df)
    
    # Save the output to a file
    output_df.to_csv(args.output_file, sep='\t', index=False)

if __name__ == "__main__":
    main() 