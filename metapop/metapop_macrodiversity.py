"""
MetaPop Macrodiversity Module

Calculates alpha and beta diversity metrics, and performs ordination analyses.
Replaces MetaPop_Macrodiversity.R.
"""

import os
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
import warnings

from .metapop_diversity_utils import (
    calculate_alpha_diversity,
    calculate_beta_diversity,
    clr_transform,
    vegdist,
    perform_pca,
    perform_pcoa,
    perform_nmds,
    specnumber,
    diversity,
    fisher_alpha,
    chao1,
    ace,
    pielou_j
)


def read_normalization_file(norm_file: str) -> pd.DataFrame:
    """
    Read and validate normalization file.

    Args:
        norm_file: Path to tab-separated normalization file

    Returns:
        DataFrame with sample names and normalization factors
    """
    norm_data = pd.read_csv(norm_file, sep='\t', header=None)

    if norm_data.shape[1] < 2:
        raise ValueError("Normalization file must have at least 2 columns")

    norm_data.columns = ['sample', 'count'] + list(norm_data.columns[2:])

    # Try to parse count as numeric
    norm_data['count'] = pd.to_numeric(norm_data['count'], errors='coerce')

    # If first row is header, re-read with header
    if pd.isna(norm_data['count'].iloc[0]):
        norm_data = pd.read_csv(norm_file, sep='\t')
        norm_data.columns = ['sample', 'count'] + list(norm_data.columns[2:])

    # Calculate normalization factor
    max_count = norm_data['count'].max()
    norm_data['factor'] = norm_data['count'] / max_count

    return norm_data


def get_fasta_contig_names(fasta_file: str) -> List[str]:
    """
    Get contig names from a FASTA file.

    Args:
        fasta_file: Path to FASTA file

    Returns:
        List of contig names
    """
    names = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                name = line[1:].strip().split()[0]
                names.append(name)
    return names


def build_abundance_matrix(cov_depth_dir: str, contigs: List[str],
                           coverage_cutoff: int = 0, min_contig_length: int = 5000,
                           norm_data: pd.DataFrame = None) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build abundance matrix from coverage/depth files.

    Args:
        cov_depth_dir: Directory containing breadth and depth files
        contigs: List of contig names
        coverage_cutoff: Minimum coverage percentage
        min_contig_length: Minimum contig length for detection
        norm_data: Normalization factors DataFrame

    Returns:
        Tuple of (raw_abundance_matrix, normalized_abundance_matrix)
    """
    import glob

    depth_files = glob.glob(os.path.join(cov_depth_dir, "*_breadth_and_depth.tsv"))

    # Initialize abundance matrix
    samples = []
    for f in depth_files:
        sample_name = os.path.basename(f).replace("_breadth_and_depth.tsv", "")
        samples.append(sample_name)

    # Create empty matrix
    abundance_data = {contig: {sample: 0.0 for sample in samples} for contig in contigs}

    # Read depth files
    for depth_file in depth_files:
        sample_name = os.path.basename(depth_file).replace("_breadth_and_depth.tsv", "")

        df = pd.read_csv(depth_file, sep='\t', header=None,
                         names=['contig', 'length', 'coverage', 'depth'])

        for _, row in df.iterrows():
            contig = row['contig']
            coverage = row['coverage']
            depth = row['depth']
            length = row['length']

            # Apply filters
            if coverage >= coverage_cutoff or (length * coverage / 100) >= min_contig_length:
                if contig in abundance_data:
                    abundance_data[contig][sample_name] = depth

    # Convert to DataFrame
    raw_matrix = pd.DataFrame(abundance_data).T
    raw_matrix.index.name = 'contig'

    # Normalize if normalization data provided
    if norm_data is not None:
        norm_factors = {row['sample']: row['factor']
                       for _, row in norm_data.iterrows()}
        normalized_matrix = raw_matrix.copy()
        for col in normalized_matrix.columns:
            if col in norm_factors:
                normalized_matrix[col] = normalized_matrix[col] / norm_factors[col]
    else:
        normalized_matrix = raw_matrix.copy()

    return raw_matrix, normalized_matrix


def run_macrodiversity(output_dir: str, ref_fasta: str, ref_genes: str,
                       norm_file: str, min_cov: int = 20, min_dep: int = 10,
                       whole_genome: bool = False, genome_cov_cutoff: int = 0,
                       min_bp_len: int = 5000, threads: int = 1) -> Dict:
    """
    Run complete macrodiversity analysis.

    Args:
        output_dir: Base output directory
        ref_fasta: Path to reference FASTA file
        ref_genes: Path to reference genes file
        norm_file: Path to normalization file
        min_cov: Minimum coverage threshold
        min_dep: Minimum depth threshold
        whole_genome: Treat references as whole genomes
        genome_cov_cutoff: Genome coverage cutoff
        min_bp_len: Minimum bases for detection
        threads: Number of threads

    Returns:
        Dictionary with all results
    """
    base_path = os.path.normpath(output_dir)
    metapop_dir = os.path.join(base_path, "MetaPop")
    output_path = os.path.join(metapop_dir, "11.Macrodiversity")
    viz_path = os.path.join(metapop_dir, "12.Visualizations")

    os.makedirs(output_path, exist_ok=True)
    os.makedirs(viz_path, exist_ok=True)

    results = {}

    # Get contig names
    contigs = get_fasta_contig_names(ref_fasta)
    print(f"Found {len(contigs)} contigs in reference")

    # Read normalization file
    norm_data = read_normalization_file(norm_file)
    print(f"Loaded normalization factors for {len(norm_data)} samples")

    # Build abundance matrix
    cov_depth_dir = os.path.join(metapop_dir, "03.Breadth_and_Depth")
    raw_matrix, norm_matrix = build_abundance_matrix(
        cov_depth_dir, contigs,
        coverage_cutoff=genome_cov_cutoff,
        min_contig_length=min_bp_len,
        norm_data=norm_data
    )

    # Save abundance tables
    raw_matrix.to_csv(os.path.join(output_path, "raw_abundances_table.tsv"), sep='\t')
    norm_matrix.to_csv(os.path.join(output_path, "normalized_abundances_table.tsv"), sep='\t')
    results['raw_abundance'] = raw_matrix
    results['normalized_abundance'] = norm_matrix

    # Convert to numpy for calculations (samples x species)
    abundance_array = norm_matrix.T.values

    if abundance_array.shape[0] == 0:
        print("No samples found. Cannot calculate diversity metrics.")
        return results

    # Calculate alpha diversity
    print("Calculating alpha diversity...")
    alpha_div = calculate_alpha_diversity(abundance_array)
    alpha_div.index = norm_matrix.columns
    alpha_div.to_csv(os.path.join(output_path, "Alpha_diversity_stats.tsv"), sep='\t')
    results['alpha_diversity'] = alpha_div
    print(f"Alpha diversity calculated for {len(alpha_div)} samples")

    # Calculate beta diversity
    print("Calculating beta diversity...")
    beta_results = calculate_beta_diversity(abundance_array,
                                            methods=['bray', 'jaccard', 'clr_euclidean'])

    for method, dist_matrix in beta_results.items():
        df = pd.DataFrame(dist_matrix, index=norm_matrix.columns, columns=norm_matrix.columns)
        df.to_csv(os.path.join(output_path, f"Beta_diversity_{method}_distances.tsv"), sep='\t')
        results[f'beta_{method}'] = df

    # Perform ordination analyses if we have enough samples
    n_samples = abundance_array.shape[0]
    if n_samples >= 3:
        print("Performing ordination analyses...")

        for method_name, dist_matrix in beta_results.items():
            try:
                # PCA
                pca_coords, pca_var = perform_pca(dist_matrix)
                pca_df = pd.DataFrame(pca_coords[:, :min(2, pca_coords.shape[1])],
                                     index=norm_matrix.columns,
                                     columns=['PC1', 'PC2'] if pca_coords.shape[1] >= 2 else ['PC1'])
                pca_df['var_explained'] = pca_var[:len(pca_df.columns)]
                results[f'pca_{method_name}'] = pca_df

                # PCoA
                pcoa_coords, pcoa_eigenvalues = perform_pcoa(dist_matrix)
                pcoa_df = pd.DataFrame(pcoa_coords[:, :min(2, pcoa_coords.shape[1])],
                                      index=norm_matrix.columns,
                                      columns=['MDS1', 'MDS2'] if pcoa_coords.shape[1] >= 2 else ['MDS1'])
                results[f'pcoa_{method_name}'] = pcoa_df

                # NMDS
                nmds_coords, stress = perform_nmds(dist_matrix)
                nmds_df = pd.DataFrame(nmds_coords,
                                      index=norm_matrix.columns,
                                      columns=['NMDS1', 'NMDS2'])
                nmds_df['stress'] = stress
                results[f'nmds_{method_name}'] = nmds_df

            except Exception as e:
                print(f"Warning: Ordination failed for {method_name}: {e}")
    else:
        print(f"Only {n_samples} samples found. Need at least 3 for ordination analyses.")

    print("Macrodiversity analysis complete.")
    return results


# =============================================================================
# Visualization Functions (matplotlib-based)
# =============================================================================

def plot_alpha_diversity(alpha_div: pd.DataFrame, output_path: str) -> None:
    """
    Create alpha diversity scatter plots.

    Args:
        alpha_div: Alpha diversity DataFrame
        output_path: Path to save plots
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')
    except ImportError:
        print("matplotlib not available. Skipping alpha diversity plots.")
        return

    metrics = ['Richness', 'Shannons_H', 'Simpson', 'InvSimpson', 'Fisher', 'Pielous_J', 'Chao1', 'ACE']

    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.flatten()

    for i, metric in enumerate(metrics):
        if metric not in alpha_div.columns:
            continue

        ax = axes[i]
        values = alpha_div[metric].dropna()

        ax.scatter(range(len(values)), values, c='gray', edgecolor='black', s=50)
        ax.axhline(values.median(), color='red', linestyle='--', label='Median')
        ax.axhline(values.mean(), color='red', linestyle='-', label='Mean')

        ax.set_xlabel('Samples')
        ax.set_ylabel(metric)
        ax.set_title(metric)

        # Rotate x-axis labels
        ax.set_xticks(range(len(values)))
        ax.set_xticklabels(values.index, rotation=90, fontsize=6)

    plt.tight_layout()
    plt.savefig(os.path.join(output_path, "alpha_diversity_scatterplots.pdf"), dpi=150)
    plt.close()


def plot_ordination(coords: pd.DataFrame, title: str, output_path: str,
                    color_by: pd.Series = None, xlabel: str = 'Axis 1',
                    ylabel: str = 'Axis 2') -> None:
    """
    Create ordination plot (PCA, PCoA, or NMDS).

    Args:
        coords: Coordinates DataFrame with first two columns as axes
        title: Plot title
        output_path: Path to save plot
        color_by: Optional series to color points by
        xlabel: X-axis label
        ylabel: Y-axis label
    """
    try:
        import matplotlib.pyplot as plt
        import matplotlib
        matplotlib.use('Agg')
    except ImportError:
        print("matplotlib not available. Skipping ordination plots.")
        return

    fig, ax = plt.subplots(figsize=(9, 9))

    x_col = coords.columns[0]
    y_col = coords.columns[1] if len(coords.columns) > 1 else coords.columns[0]

    if color_by is not None:
        scatter = ax.scatter(coords[x_col], coords[y_col], c=color_by,
                           cmap='terrain', s=100, edgecolor='black')
        plt.colorbar(scatter, label='Species Richness')
    else:
        ax.scatter(coords[x_col], coords[y_col], c='gray', s=100, edgecolor='black')

    # Add labels
    for idx, row in coords.iterrows():
        ax.annotate(idx, (row[x_col], row[y_col]), fontsize=8,
                   xytext=(5, 5), textcoords='offset points')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_aspect('equal')

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()


def plot_heatmap(matrix: pd.DataFrame, output_path: str, title: str = 'Abundance Heatmap') -> None:
    """
    Create clustered heatmap of abundance matrix.

    Args:
        matrix: Abundance matrix (contigs x samples)
        output_path: Path to save plot
        title: Plot title
    """
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        import matplotlib
        matplotlib.use('Agg')
    except ImportError:
        print("matplotlib/seaborn not available. Skipping heatmap.")
        return

    # Remove all-zero rows
    matrix_nozero = matrix[matrix.sum(axis=1) > 0]

    if len(matrix_nozero) == 0:
        print("No non-zero rows in matrix. Skipping heatmap.")
        return

    # Limit to 75th percentile for better visualization
    max_val = np.percentile(matrix_nozero.values[matrix_nozero.values > 0], 75)
    matrix_clipped = matrix_nozero.clip(upper=max_val)

    try:
        plt.figure(figsize=(12, 10))
        sns.clustermap(matrix_clipped, cmap='Spectral_r', yticklabels=False,
                      xticklabels=True, figsize=(12, 10))
        plt.savefig(output_path, dpi=150)
        plt.close()
    except Exception as e:
        print(f"Heatmap failed: {e}")


def generate_macrodiversity_visualizations(results: Dict, output_path: str) -> None:
    """
    Generate all macrodiversity visualizations.

    Args:
        results: Dictionary from run_macrodiversity()
        output_path: Path to save visualizations
    """
    os.makedirs(output_path, exist_ok=True)

    # Alpha diversity plots
    if 'alpha_diversity' in results:
        plot_alpha_diversity(results['alpha_diversity'], output_path)

    # Get richness for coloring
    richness = None
    if 'alpha_diversity' in results and 'Richness' in results['alpha_diversity'].columns:
        richness = results['alpha_diversity']['Richness']

    # Ordination plots
    for key in results:
        if key.startswith('pca_'):
            method = key.replace('pca_', '')
            coords = results[key]
            var_exp = coords.get('var_explained', [0, 0])
            plot_ordination(
                coords[['PC1', 'PC2']] if 'PC2' in coords.columns else coords[['PC1']],
                f'PCA - {method}',
                os.path.join(output_path, f'PCA_{method}_plot.pdf'),
                color_by=richness,
                xlabel=f'PC1 ({var_exp[0]*100:.1f}%)' if len(var_exp) > 0 else 'PC1',
                ylabel=f'PC2 ({var_exp[1]*100:.1f}%)' if len(var_exp) > 1 else 'PC2'
            )

        elif key.startswith('pcoa_'):
            method = key.replace('pcoa_', '')
            coords = results[key]
            plot_ordination(
                coords[['MDS1', 'MDS2']] if 'MDS2' in coords.columns else coords[['MDS1']],
                f'PCoA - {method}',
                os.path.join(output_path, f'PCoA_{method}_plot.pdf'),
                color_by=richness,
                xlabel='PCo1',
                ylabel='PCo2'
            )

        elif key.startswith('nmds_'):
            method = key.replace('nmds_', '')
            coords = results[key]
            stress = coords['stress'].iloc[0] if 'stress' in coords.columns else 0
            plot_ordination(
                coords[['NMDS1', 'NMDS2']],
                f'NMDS - {method} (stress={stress:.3f})',
                os.path.join(output_path, f'NMDS_{method}_plot.pdf'),
                color_by=richness,
                xlabel='NMDS1',
                ylabel='NMDS2'
            )

    # Heatmap
    if 'normalized_abundance' in results:
        plot_heatmap(
            results['normalized_abundance'],
            os.path.join(output_path, 'normalized_abundances_heatmap.pdf')
        )

    print("Macrodiversity visualizations complete.")
