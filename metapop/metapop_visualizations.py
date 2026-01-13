"""
MetaPop Visualization Module

Unified visualization module combining preprocessing summaries, microdiversity,
and codon bias visualizations. Replaces R scripts:
- MetaPop_Preprocessing_Summaries.R
- MetaPop_Microdiversity_Visualizations.R
- MetaPop_Codon_Bias_Viz.R
"""

import os
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
import glob
import warnings

from .metapop_plotting_utils import (
    MATPLOTLIB_AVAILABLE, SEABORN_AVAILABLE,
    get_brewer_palette, get_pnps_colormap,
    draw_gene_arrow, plot_gene_track,
    plot_depth_coverage, plot_clustered_heatmap,
    plot_fst_heatmap, plot_scatter_with_stats,
    plot_donut_chart, plot_stacked_bar,
    create_figure_grid, save_figure
)

if MATPLOTLIB_AVAILABLE:
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('Agg')

if SEABORN_AVAILABLE:
    import seaborn as sns


# =============================================================================
# Preprocessing Visualizations
# =============================================================================

def visualize_preprocessing(output_dir: str, original_bams_dir: str,
                           sample_names: List[str]) -> None:
    """
    Generate preprocessing summary visualizations.

    Args:
        output_dir: Base output directory
        original_bams_dir: Directory with original BAM files
        sample_names: List of sample names
    """
    if not MATPLOTLIB_AVAILABLE:
        print("matplotlib not available. Skipping preprocessing visualizations.")
        return

    base_path = os.path.normpath(output_dir)
    metapop_dir = os.path.join(base_path, "MetaPop")
    viz_path = os.path.join(metapop_dir, "12.Visualizations")
    os.makedirs(viz_path, exist_ok=True)

    # Read breadth/depth data
    bd_path = os.path.join(metapop_dir, "03.Breadth_and_Depth")
    bd_files = glob.glob(os.path.join(bd_path, "*_breadth_and_depth.tsv"))

    if not bd_files:
        print("No breadth/depth files found. Skipping preprocessing viz.")
        return

    # Aggregate statistics
    stats = []
    for f in bd_files:
        sample = os.path.basename(f).replace("_breadth_and_depth.tsv", "")
        df = pd.read_csv(f, sep='\t', header=None,
                        names=['contig', 'length', 'coverage', 'depth'])
        stats.append({
            'sample': sample,
            'num_contigs': len(df),
            'total_length': df['length'].sum(),
            'mean_coverage': df['coverage'].mean(),
            'mean_depth': df['depth'].mean(),
            'contigs_passing': (df['coverage'] >= 20).sum()
        })

    stats_df = pd.DataFrame(stats)

    # Create summary figure
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Coverage distribution
    ax = axes[0, 0]
    for f in bd_files[:5]:  # Limit to first 5 for readability
        sample = os.path.basename(f).replace("_breadth_and_depth.tsv", "")
        df = pd.read_csv(f, sep='\t', header=None,
                        names=['contig', 'length', 'coverage', 'depth'])
        ax.hist(df['coverage'], bins=20, alpha=0.5, label=sample)
    ax.set_xlabel('Coverage (%)')
    ax.set_ylabel('Count')
    ax.set_title('Coverage Distribution')
    ax.legend(fontsize=6)

    # Depth distribution
    ax = axes[0, 1]
    for f in bd_files[:5]:
        sample = os.path.basename(f).replace("_breadth_and_depth.tsv", "")
        df = pd.read_csv(f, sep='\t', header=None,
                        names=['contig', 'length', 'coverage', 'depth'])
        ax.hist(df['depth'], bins=20, alpha=0.5, label=sample)
    ax.set_xlabel('Depth')
    ax.set_ylabel('Count')
    ax.set_title('Depth Distribution')
    ax.legend(fontsize=6)

    # Contigs passing per sample
    ax = axes[1, 0]
    ax.bar(stats_df['sample'], stats_df['contigs_passing'], color='steelblue', edgecolor='black')
    ax.set_xlabel('Sample')
    ax.set_ylabel('Contigs Passing Filters')
    ax.set_title('Contigs Passing per Sample')
    ax.tick_params(axis='x', rotation=45)

    # Mean depth per sample
    ax = axes[1, 1]
    ax.bar(stats_df['sample'], stats_df['mean_depth'], color='coral', edgecolor='black')
    ax.set_xlabel('Sample')
    ax.set_ylabel('Mean Depth')
    ax.set_title('Mean Depth per Sample')
    ax.tick_params(axis='x', rotation=45)

    plt.tight_layout()
    plt.savefig(os.path.join(viz_path, 'preprocessing_summary.pdf'), dpi=150)
    plt.close()

    print("Preprocessing visualizations complete.")


# =============================================================================
# Microdiversity Visualizations
# =============================================================================

def visualize_microdiversity(output_dir: str, plot_all: bool = False,
                            snp_scale: str = 'local') -> None:
    """
    Generate microdiversity visualizations.

    Args:
        output_dir: Base output directory
        plot_all: Plot all genomes (not just top 20)
        snp_scale: 'local', 'global', or 'both'
    """
    if not MATPLOTLIB_AVAILABLE:
        print("matplotlib not available. Skipping microdiversity visualizations.")
        return

    base_path = os.path.normpath(output_dir)
    metapop_dir = os.path.join(base_path, "MetaPop")
    micro_path = os.path.join(metapop_dir, "10.Microdiversity")
    viz_path = os.path.join(metapop_dir, "12.Visualizations")
    os.makedirs(viz_path, exist_ok=True)

    # Check for required files
    codon_pos_file = os.path.join(micro_path, "global_codon_position_summary.tsv")
    gene_micro_file = os.path.join(micro_path, "global_gene_microdiversity.tsv")

    # Codon position summary
    if os.path.exists(codon_pos_file):
        codon_pos = pd.read_csv(codon_pos_file, sep='\t')

        fig, ax = plt.subplots(figsize=(8, 6))
        totals = codon_pos[['first_pos', 'second_pos', 'third_pos']].sum()
        bars = ax.bar(['1st Position', '2nd Position', '3rd Position'],
                     totals.values, color=['#e41a1c', '#377eb8', '#4daf4a'],
                     edgecolor='black')
        ax.set_ylabel('SNP Count')
        ax.set_title('SNPs by Codon Position')

        # Add warning if 3rd position is not highest
        if totals['third_pos'] < max(totals['first_pos'], totals['second_pos']):
            ax.text(0.5, 0.95, 'Warning: Unusual codon position distribution',
                   transform=ax.transAxes, ha='center', color='red', fontsize=10)

        plt.tight_layout()
        plt.savefig(os.path.join(viz_path, 'codon_position_summary.pdf'), dpi=150)
        plt.close()

    # Gene-level microdiversity
    if os.path.exists(gene_micro_file):
        gene_micro = pd.read_csv(gene_micro_file, sep='\t')

        if len(gene_micro) > 0:
            # pN/pS distribution
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))

            # pN/pS histogram
            ax = axes[0, 0]
            pnps_values = gene_micro['pNpS_ratio'].dropna()
            pnps_values = pnps_values[~np.isinf(pnps_values)]
            if len(pnps_values) > 0:
                ax.hist(pnps_values, bins=50, color='steelblue', edgecolor='black')
                ax.axvline(1.0, color='red', linestyle='--', label='Neutral (pN/pS=1)')
                ax.set_xlabel('pN/pS Ratio')
                ax.set_ylabel('Count')
                ax.set_title('pN/pS Distribution')
                ax.legend()

            # Pi vs Theta scatter
            ax = axes[0, 1]
            pi_vals = gene_micro['pi'].dropna()
            theta_vals = gene_micro['theta'].dropna()
            if len(pi_vals) > 0 and len(theta_vals) > 0:
                ax.scatter(theta_vals[:len(pi_vals)], pi_vals[:len(theta_vals)],
                          alpha=0.5, c='steelblue')
                ax.plot([0, max(theta_vals)], [0, max(theta_vals)], 'r--', label='Pi = Theta')
                ax.set_xlabel('Theta (Watterson)')
                ax.set_ylabel('Pi (Nucleotide Diversity)')
                ax.set_title('Pi vs Theta')
                ax.legend()

            # Tajima's D distribution
            ax = axes[1, 0]
            taj_d = gene_micro['taj_D'].dropna()
            taj_d = taj_d[~np.isinf(taj_d)]
            if len(taj_d) > 0:
                ax.hist(taj_d, bins=50, color='coral', edgecolor='black')
                ax.axvline(0, color='black', linestyle='--')
                ax.axvline(-2, color='red', linestyle=':', label='Purifying selection')
                ax.axvline(2, color='blue', linestyle=':', label='Balancing selection')
                ax.set_xlabel("Tajima's D")
                ax.set_ylabel('Count')
                ax.set_title("Tajima's D Distribution")
                ax.legend()

            # SNPs per gene
            ax = axes[1, 1]
            snp_counts = gene_micro['num_snps']
            ax.hist(snp_counts, bins=50, color='green', edgecolor='black')
            ax.set_xlabel('Number of SNPs')
            ax.set_ylabel('Gene Count')
            ax.set_title('SNPs per Gene')

            plt.tight_layout()
            plt.savefig(os.path.join(viz_path, 'microdiversity_summary.pdf'), dpi=150)
            plt.close()

    # FST heatmaps
    fst_file = os.path.join(micro_path, "fst_results.tsv")
    if os.path.exists(fst_file):
        fst_data = pd.read_csv(fst_file, sep='\t', index_col=0)
        plot_fst_heatmap(fst_data, os.path.join(viz_path, 'fst_heatmap.pdf'))

    print("Microdiversity visualizations complete.")


# =============================================================================
# Codon Bias Visualizations
# =============================================================================

def visualize_codon_bias(output_dir: str) -> None:
    """
    Generate codon bias visualizations.

    Args:
        output_dir: Base output directory
    """
    if not MATPLOTLIB_AVAILABLE:
        print("matplotlib not available. Skipping codon bias visualizations.")
        return

    base_path = os.path.normpath(output_dir)
    metapop_dir = os.path.join(base_path, "MetaPop")
    cb_path = os.path.join(metapop_dir, "08.Codon_Bias")
    viz_path = os.path.join(metapop_dir, "12.Visualizations")
    os.makedirs(viz_path, exist_ok=True)

    # Gene Euclidean distances
    euc_dist_file = os.path.join(cb_path, "gene_euclidean_distances.tsv")
    if os.path.exists(euc_dist_file):
        euc_dist = pd.read_csv(euc_dist_file, sep='\t')

        fig, axes = plt.subplots(1, 2, figsize=(12, 5))

        # Euclidean distance distribution
        ax = axes[0]
        ax.hist(euc_dist['euc_dist'], bins=50, color='purple', edgecolor='black', alpha=0.7)
        ax.set_xlabel('Euclidean Distance from Mean')
        ax.set_ylabel('Gene Count')
        ax.set_title('Codon Bias Distance Distribution')

        # Outlier vs non-outlier
        ax = axes[1]
        if 'outlier_status' in euc_dist.columns:
            outlier_counts = euc_dist['outlier_status'].value_counts()
            colors = ['#2ecc71', '#e74c3c']
            ax.pie(outlier_counts.values, labels=outlier_counts.index,
                  colors=colors, autopct='%1.1f%%',
                  wedgeprops=dict(edgecolor='black'))
            ax.set_title('Outlier Gene Distribution')

        plt.tight_layout()
        plt.savefig(os.path.join(viz_path, 'codon_bias_summary.pdf'), dpi=150)
        plt.close()

    # Per-contig codon bias
    iqr_file = os.path.join(cb_path, "gene_IQR_and_mean.tsv")
    if os.path.exists(iqr_file):
        iqr_data = pd.read_csv(iqr_file, sep='\t')

        if len(iqr_data) > 0:
            fig, ax = plt.subplots(figsize=(10, 6))

            # Polar plot for top contigs
            n_contigs = min(20, len(iqr_data))
            top_contigs = iqr_data.nlargest(n_contigs, 'mu')

            angles = np.linspace(0, 2 * np.pi, n_contigs, endpoint=False)
            radii = top_contigs['mu'].values

            bars = ax.bar(range(n_contigs), radii, color='steelblue', edgecolor='black')
            ax.set_xticks(range(n_contigs))
            ax.set_xticklabels(top_contigs['parent_contig'].values, rotation=45, ha='right', fontsize=6)
            ax.set_ylabel('Mean Euclidean Distance')
            ax.set_title('Codon Bias by Contig (Top 20)')

            plt.tight_layout()
            plt.savefig(os.path.join(viz_path, 'codon_bias_by_contig.pdf'), dpi=150)
            plt.close()

    print("Codon bias visualizations complete.")


# =============================================================================
# Main Visualization Entry Point
# =============================================================================

def run_all_visualizations(output_dir: str, original_bams_dir: str = None,
                          sample_names: List[str] = None,
                          plot_all: bool = False, snp_scale: str = 'local') -> None:
    """
    Run all visualization routines.

    Args:
        output_dir: Base output directory
        original_bams_dir: Directory with original BAM files
        sample_names: List of sample names
        plot_all: Plot all genomes for microdiversity
        snp_scale: SNP scale for microdiversity viz
    """
    print("Generating visualizations...")

    # Preprocessing
    if original_bams_dir and sample_names:
        try:
            visualize_preprocessing(output_dir, original_bams_dir, sample_names)
        except Exception as e:
            print(f"Preprocessing visualization failed: {e}")

    # Microdiversity
    try:
        visualize_microdiversity(output_dir, plot_all, snp_scale)
    except Exception as e:
        print(f"Microdiversity visualization failed: {e}")

    # Codon bias
    try:
        visualize_codon_bias(output_dir)
    except Exception as e:
        print(f"Codon bias visualization failed: {e}")

    print("All visualizations complete.")
