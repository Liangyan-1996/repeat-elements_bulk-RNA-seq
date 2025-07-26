import os
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.utils import load_example_data

def load_count_data(count_dir):
    """
    Load count data from .cntTable files in the specified directory.
    File names are expected to follow the pattern: {group}_rep{number}.cntTable
    """
    count_files = [f for f in os.listdir(count_dir) if f.endswith('.cntTable')]
    
    # Dictionary to store counts for each group
    group_counts = {}
    sample_info = []
    
    for file in count_files:
        # Extract group name from file name
        group = file.split('_rep')[0]
        
        # Read count data
        file_path = os.path.join(count_dir, file)
        counts = pd.read_csv(file_path, sep='\t', header=0)
        counts.columns = ['gene_id', 'counts']
        
        # Add to group counts
        if group not in group_counts:
            group_counts[group] = counts[['gene_id']]
        
        # Add counts as a new column
        counts_col = counts['counts']
        counts_col.name = file.replace('.cntTable', '')
        group_counts[group] = pd.concat([group_counts[group], counts_col], axis=1)
        
        # Add sample info
        sample_info.append({'sample': file.replace('.cntTable', ''), 'group': group})
    
    # Combine all count data into a single matrix
    all_counts = pd.concat(group_counts.values(), axis=1)
    # Remove duplicate gene_id columns
    all_counts = all_counts.loc[:,~all_counts.columns.duplicated()]
    
    # Create sample info dataframe
    sample_df = pd.DataFrame(sample_info)
    sample_df.set_index('sample', inplace=True)
    
    return all_counts, sample_df

def run_differential_expression(counts_df, sample_df, control_group='Nontargeting_control'):
    """
    Run differential expression analysis using pydeseq2
    """
    # Get all groups
    groups = sample_df['group'].unique()
    
    # Run analysis for each group vs control
    results = {}
    for group in groups:
        if group != control_group:
            # Create contrast
            contrast_samples = sample_df[sample_df['group'].isin([control_group, group])]
            contrast_counts = counts_df[contrast_samples.index]
            
            # Create design matrix
            design_matrix = contrast_samples.copy()
            
            # Initialize DeseqDataSet
            dds = DeseqDataSet(
                counts=contrast_counts,
                metadata=design_matrix,
                design_factors="group",
                refit_cooks=True,
                n_cpus=16
            )
            
            # Run differential expression analysis
            dds.deseq2()
            
            # Get results
            stat_res = DeseqStats(dds, contrast=["group", group, control_group])
            stat_res.summary()
            stat_res.lfc_shrink(coeff="condition[T.B]")
            results[group] = stat_res.results_df
    
    return results

def save_results(results, output_dir='DESeq2_results'):
    """
    Save differential expression results to CSV files
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for group, result_df in results.items():
        output_file = os.path.join(output_dir, f"{group}_vs_Nontargeting_control.csv")
        result_df.to_csv(output_file)
        print(f"Results for {group} saved to {output_file}")

def main():
    # Load count data
    counts_df, sample_df = load_count_data('TEcount')
    
    # Run differential expression analysis
    results = run_differential_expression(counts_df, sample_df)
    
    # Save results
    save_results(results)

if __name__ == "__main__":
    main()