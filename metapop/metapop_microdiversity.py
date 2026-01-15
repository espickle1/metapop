"""
MetaPop Microdiversity Module

Calculates microdiversity metrics including Tajima's D, pi (nucleotide diversity),
theta (Watterson's estimator), and pN/pS ratios. Replaces MetaPop_Microdiversity.R.
"""

import os
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional
from concurrent.futures import ProcessPoolExecutor
import warnings

from .metapop_diversity_utils import tajima_d_parameters, calculate_tajimas_d


# =============================================================================
# Codon and Amino Acid Constants
# =============================================================================

CODONS = [
    "AAA", "AAT", "AAC", "AAG", "ATA", "ATT", "ATC", "ATG",
    "ACA", "ACT", "ACC", "ACG", "AGA", "AGT", "AGC", "AGG",
    "TAA", "TAT", "TAC", "TAG", "TTA", "TTT", "TTC", "TTG",
    "TCA", "TCT", "TCC", "TCG", "TGA", "TGT", "TGC", "TGG",
    "CAA", "CAT", "CAC", "CAG", "CTA", "CTT", "CTC", "CTG",
    "CCA", "CCT", "CCC", "CCG", "CGA", "CGT", "CGC", "CGG",
    "GAA", "GAT", "GAC", "GAG", "GTA", "GTT", "GTC", "GTG",
    "GCA", "GCT", "GCC", "GCG", "GGA", "GGT", "GGC", "GGG"
]

CODON_TO_AA = {
    "AAA": "K", "AAT": "N", "AAC": "N", "AAG": "K",
    "ATA": "I", "ATT": "I", "ATC": "I", "ATG": "M",
    "ACA": "T", "ACT": "T", "ACC": "T", "ACG": "T",
    "AGA": "R", "AGT": "S", "AGC": "S", "AGG": "R",
    "TAA": "STOP", "TAT": "Y", "TAC": "Y", "TAG": "STOP",
    "TTA": "L", "TTT": "F", "TTC": "F", "TTG": "L",
    "TCA": "S", "TCT": "S", "TCC": "S", "TCG": "S",
    "TGA": "STOP", "TGT": "C", "TGC": "C", "TGG": "W",
    "CAA": "Q", "CAT": "H", "CAC": "H", "CAG": "Q",
    "CTA": "L", "CTT": "L", "CTC": "L", "CTG": "L",
    "CCA": "P", "CCT": "P", "CCC": "P", "CCG": "P",
    "CGA": "R", "CGT": "R", "CGC": "R", "CGG": "R",
    "GAA": "E", "GAT": "D", "GAC": "D", "GAG": "E",
    "GTA": "V", "GTT": "V", "GTC": "V", "GTG": "V",
    "GCA": "A", "GCT": "A", "GCC": "A", "GCG": "A",
    "GGA": "G", "GGT": "G", "GGC": "G", "GGG": "G"
}

# Expected N (nonsynonymous) values per codon
EXP_N_CONSTRUCT = [
    8/3, 8/3, 8/3, 8/3, 7/3, 7/3, 7/3, 3,
    2, 2, 2, 2, 7/3, 8/3, 8/3, 7/3,
    7/3, 8/3, 8/3, 8/3, 7/3, 8/3, 8/3, 7/3,
    2, 2, 2, 2, 8/3, 8/3, 8/3, 3,
    8/3, 8/3, 8/3, 8/3, 5/3, 1, 1, 5/3,
    2, 2, 2, 2, 5/3, 2, 2, 5/3,
    8/3, 8/3, 8/3, 8/3, 2, 2, 2, 2,
    2, 2, 2, 2, 2, 2, 2, 2
]

# Expected S (synonymous) values per codon
EXP_S_CONSTRUCT = [3 - n for n in EXP_N_CONSTRUCT]


def get_codon_to_index() -> Dict[str, int]:
    """Get mapping from codon to index."""
    return {codon: i for i, codon in enumerate(CODONS)}


def get_reverse_complement(codon: str) -> str:
    """Get reverse complement of a codon."""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement.get(base, 'N') for base in reversed(codon))


# =============================================================================
# Precompute Tajima's D Parameters
# =============================================================================

def precompute_tajima_params(max_sample_size: int = 50) -> Dict[int, Dict[str, float]]:
    """
    Precompute Tajima's D parameters for sample sizes 1 to max_sample_size.

    Args:
        max_sample_size: Maximum sample size to precompute

    Returns:
        Dictionary mapping sample size to parameter dict
    """
    return {n: tajima_d_parameters(n) for n in range(1, max_sample_size + 1)}


# =============================================================================
# Pi (Nucleotide Diversity) Calculation
# =============================================================================

def calculate_pi_per_site(a_ct: int, t_ct: int, c_ct: int, g_ct: int) -> float:
    """
    Calculate nucleotide diversity (pi) for a single site.

    Pi = 2 * sum(pairwise differences) / (n * (n-1))

    Args:
        a_ct, t_ct, c_ct, g_ct: Counts of each nucleotide

    Returns:
        Pi value for this site
    """
    depth = a_ct + t_ct + c_ct + g_ct

    if depth <= 1:
        return 0.0

    # Number of pairwise differences
    # For each pair of different nucleotides, count a_ct * other_ct
    pairwise_diff = (
        a_ct * (t_ct + c_ct + g_ct) +
        t_ct * (c_ct + g_ct) +
        c_ct * g_ct
    )

    # Total number of pairs
    total_pairs = depth * (depth - 1) / 2

    return pairwise_diff / total_pairs if total_pairs > 0 else 0.0


def calculate_pi_vectorized(counts: pd.DataFrame) -> np.ndarray:
    """
    Calculate pi for multiple sites vectorized.

    Args:
        counts: DataFrame with columns a_ct, t_ct, c_ct, g_ct

    Returns:
        Array of pi values
    """
    a = counts['a_ct'].values if 'a_ct' in counts else counts['sub_samp_a'].values
    t = counts['t_ct'].values if 't_ct' in counts else counts['sub_samp_t'].values
    c = counts['c_ct'].values if 'c_ct' in counts else counts['sub_samp_c'].values
    g = counts['g_ct'].values if 'g_ct' in counts else counts['sub_samp_g'].values

    depth = a + t + c + g

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pairwise_diff = (
            a * (t + c + g) +
            t * (c + g) +
            c * g
        )
        total_pairs = (depth * (depth - 1)) / 2
        pi = np.where(total_pairs > 0, pairwise_diff / total_pairs, 0.0)

    return pi


# =============================================================================
# Subsampling Functions
# =============================================================================

def subsample_counts(row: pd.Series, sub_samp: int) -> Dict[str, int]:
    """
    Subsample nucleotide counts to normalize depth.

    Args:
        row: Series with a_ct, t_ct, c_ct, g_ct columns
        sub_samp: Target subsample size

    Returns:
        Dictionary with subsampled counts
    """
    depth = row['a_ct'] + row['t_ct'] + row['c_ct'] + row['g_ct']

    if depth <= sub_samp:
        return {
            'sub_samp_a': row['a_ct'],
            'sub_samp_t': row['t_ct'],
            'sub_samp_c': row['c_ct'],
            'sub_samp_g': row['g_ct'],
            'sub_samp_depth': depth
        }

    prop = sub_samp / depth
    sub_a = round(row['a_ct'] * prop)
    sub_t = round(row['t_ct'] * prop)
    sub_c = round(row['c_ct'] * prop)
    sub_g = round(row['g_ct'] * prop)

    return {
        'sub_samp_a': sub_a,
        'sub_samp_t': sub_t,
        'sub_samp_c': sub_c,
        'sub_samp_g': sub_g,
        'sub_samp_depth': sub_a + sub_t + sub_c + sub_g
    }


# =============================================================================
# pN/pS Calculation Functions
# =============================================================================

def classify_snp_ns(original_codon: str, snp_base: str, pos_in_codon: int) -> str:
    """
    Classify a SNP as synonymous (S) or nonsynonymous (N).

    Args:
        original_codon: The reference codon (3 bases)
        snp_base: The alternate base
        pos_in_codon: Position in codon (1, 2, or 3)

    Returns:
        'N' for nonsynonymous, 'S' for synonymous
    """
    if original_codon not in CODON_TO_AA:
        return 'N'  # Default to nonsynonymous for invalid codons

    original_aa = CODON_TO_AA[original_codon]

    # Create mutant codon
    codon_list = list(original_codon)
    codon_list[pos_in_codon - 1] = snp_base
    mutant_codon = ''.join(codon_list)

    if mutant_codon not in CODON_TO_AA:
        return 'N'

    mutant_aa = CODON_TO_AA[mutant_codon]

    return 'S' if original_aa == mutant_aa else 'N'


def calculate_gene_exp_ns(gene_sequence: str) -> Tuple[float, float]:
    """
    Calculate expected N and S for a gene sequence.

    Args:
        gene_sequence: DNA sequence of the gene

    Returns:
        Tuple of (expected_N, expected_S)
    """
    codon_to_idx = get_codon_to_index()
    exp_n = 0.0
    exp_s = 0.0

    # Split sequence into codons
    for i in range(0, len(gene_sequence) - 2, 3):
        codon = gene_sequence[i:i+3].upper()
        if codon in codon_to_idx:
            idx = codon_to_idx[codon]
            exp_n += EXP_N_CONSTRUCT[idx]
            exp_s += EXP_S_CONSTRUCT[idx]

    return exp_n, exp_s


def calculate_pnps_ratio(obs_n: int, obs_s: int, exp_n: float, exp_s: float) -> float:
    """
    Calculate pN/pS ratio.

    pN/pS = (obs_N / obs_S) / (exp_N / exp_S)

    Args:
        obs_n: Observed nonsynonymous SNPs
        obs_s: Observed synonymous SNPs
        exp_n: Expected nonsynonymous sites
        exp_s: Expected synonymous sites

    Returns:
        pN/pS ratio
    """
    if obs_s == 0 or exp_s == 0 or exp_n == 0:
        return np.nan

    return (obs_n / obs_s) / (exp_n / exp_s)


# =============================================================================
# Main Microdiversity Processing Functions
# =============================================================================

def geometric_mean(x: np.ndarray) -> float:
    """Calculate geometric mean, ceiling the result."""
    x = x[x > 0]
    if len(x) == 0:
        return 0
    return np.ceil(np.exp(np.sum(np.log(x)) / len(x)))


def process_genic_snps(genic_snps_file: str, sub_samp: int = 10,
                       output_dir: str = None) -> pd.DataFrame:
    """
    Process genic SNP data and calculate microdiversity metrics.

    Args:
        genic_snps_file: Path to genic SNPs TSV file
        sub_samp: Subsample size for normalization
        output_dir: Output directory for results

    Returns:
        DataFrame with processed microdiversity data
    """
    # Read genic SNPs
    microdiv_data = pd.read_csv(genic_snps_file, sep='\t')

    # Calculate depth
    microdiv_data['depth'] = (
        microdiv_data['a_ct'] + microdiv_data['t_ct'] +
        microdiv_data['c_ct'] + microdiv_data['g_ct']
    )

    # Filter out positions with depth <= 1
    microdiv_data = microdiv_data[microdiv_data['depth'] > 1].copy()

    if len(microdiv_data) == 0:
        return pd.DataFrame()

    # Subsample counts
    microdiv_data['sample_prop'] = sub_samp / microdiv_data['depth']

    for col, orig in [('sub_samp_a', 'a_ct'), ('sub_samp_t', 't_ct'),
                      ('sub_samp_c', 'c_ct'), ('sub_samp_g', 'g_ct')]:
        microdiv_data[col] = np.where(
            microdiv_data['sample_prop'] < 1,
            np.round(microdiv_data[orig] * microdiv_data['sample_prop']),
            microdiv_data[orig]
        ).astype(int)

    microdiv_data['sub_samp_depth'] = (
        microdiv_data['sub_samp_a'] + microdiv_data['sub_samp_t'] +
        microdiv_data['sub_samp_c'] + microdiv_data['sub_samp_g']
    )

    # Calculate pi for each site
    microdiv_data['pi'] = calculate_pi_vectorized(microdiv_data)

    return microdiv_data


def calculate_gene_level_microdiversity(microdiv_data: pd.DataFrame,
                                        gene_assembly: pd.DataFrame,
                                        tajima_params: Dict[int, Dict],
                                        sub_samp: int = 10) -> pd.DataFrame:
    """
    Calculate gene-level microdiversity statistics.

    Args:
        microdiv_data: Processed SNP data from process_genic_snps()
        gene_assembly: Gene information with sequences
        tajima_params: Precomputed Tajima's D parameters
        sub_samp: Subsample size

    Returns:
        DataFrame with gene-level microdiversity
    """
    if len(microdiv_data) == 0:
        return pd.DataFrame()

    # Group by gene and source
    grouped = microdiv_data.groupby(['contig_gene', 'source'])

    results = []
    for (gene, source), group in grouped:
        contig = group['contig'].iloc[0]
        pi_sum = group['pi'].sum()
        num_snps = len(group)
        harmonic_num = int(geometric_mean(group['sub_samp_depth'].values))

        if harmonic_num < 3:
            continue

        # Get gene info
        gene_info = gene_assembly[gene_assembly['gene'] == gene]
        if len(gene_info) == 0:
            continue

        gene_len = gene_info['length'].iloc[0]
        exp_n = gene_info['expN'].iloc[0] if 'expN' in gene_info else 0
        exp_s = gene_info['expS'].iloc[0] if 'expS' in gene_info else 0

        # Get Tajima's D parameters
        params = tajima_params.get(harmonic_num, tajima_params[sub_samp])

        # Calculate theta (Watterson's estimator)
        theta = (num_snps / params['harmonic']) / gene_len if gene_len > 0 else 0

        # Normalize pi by gene length
        pi_normalized = pi_sum / gene_len if gene_len > 0 else 0

        # Calculate Tajima's D
        taj_d = calculate_tajimas_d(pi_sum, num_snps, params)

        # Calculate observed N/S
        obs_n = group['obsN'].sum() if 'obsN' in group else 0
        obs_s = group['obsS'].sum() if 'obsS' in group else 0

        # Calculate pN/pS
        pnps = calculate_pnps_ratio(obs_n, obs_s, exp_n, exp_s)

        results.append({
            'contig_gene': gene,
            'source': source,
            'contig': contig,
            'pi': pi_normalized,
            'num_snps': num_snps,
            'gene_len': gene_len,
            'theta': theta,
            'taj_D': taj_d,
            'expN': exp_n,
            'expS': exp_s,
            'obsN': obs_n,
            'obsS': obs_s,
            'pNpS_ratio': pnps,
            'snps_present': True
        })

    return pd.DataFrame(results)


def calculate_contig_level_microdiversity(microdiv_data: pd.DataFrame,
                                          contig_lengths: Dict[str, int],
                                          tajima_params: Dict[int, Dict],
                                          sub_samp: int = 10) -> pd.DataFrame:
    """
    Calculate contig-level microdiversity statistics.

    Args:
        microdiv_data: Processed SNP data
        contig_lengths: Dictionary of contig names to lengths
        tajima_params: Precomputed Tajima's D parameters
        sub_samp: Subsample size

    Returns:
        DataFrame with contig-level microdiversity
    """
    if len(microdiv_data) == 0:
        return pd.DataFrame()

    # Group by contig and source
    grouped = microdiv_data.groupby(['contig', 'source'])

    results = []
    for (contig, source), group in grouped:
        pi_sum = group['pi'].sum()
        num_snps = len(group)
        harmonic_num = int(geometric_mean(group['sub_samp_depth'].values))

        if harmonic_num < 3:
            continue

        contig_len = contig_lengths.get(contig, 0)
        if contig_len == 0:
            continue

        params = tajima_params.get(harmonic_num, tajima_params[sub_samp])

        theta = (num_snps / params['harmonic']) / contig_len
        pi_normalized = pi_sum / contig_len

        results.append({
            'contig': contig,
            'source': source,
            'pi': pi_normalized,
            'num_snps': num_snps,
            'contig_len': contig_len,
            'theta': theta,
            'snps_present': True
        })

    return pd.DataFrame(results)


def get_fasta_lengths(fasta_file: str) -> Dict[str, int]:
    """
    Get sequence names and lengths from a FASTA file.

    Args:
        fasta_file: Path to FASTA file

    Returns:
        Dictionary mapping sequence names to lengths
    """
    lengths = {}
    current_name = None
    current_length = 0

    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_name:
                    lengths[current_name] = current_length
                current_name = line[1:].split()[0]
                current_length = 0
            else:
                current_length += len(line)

        if current_name:
            lengths[current_name] = current_length

    return lengths


def calculate_codon_position_summary(microdiv_data: pd.DataFrame) -> pd.DataFrame:
    """
    Summarize SNP counts by codon position.

    Args:
        microdiv_data: Processed SNP data with pos_in_codon column

    Returns:
        DataFrame with position counts per source/contig
    """
    if 'pos_in_codon' not in microdiv_data.columns:
        return pd.DataFrame()

    grouped = microdiv_data.groupby(['source', 'contig'])

    results = []
    for (source, contig), group in grouped:
        first_pos = (group['pos_in_codon'] == 1).sum()
        second_pos = (group['pos_in_codon'] == 2).sum()
        third_pos = (group['pos_in_codon'] == 3).sum()

        results.append({
            'source': source,
            'contig': contig,
            'first_pos': first_pos,
            'second_pos': second_pos,
            'third_pos': third_pos
        })

    return pd.DataFrame(results)


# =============================================================================
# Main Entry Point
# =============================================================================

def run_microdiversity(output_dir: str, ref_fasta: str, ref_genes: str,
                       min_cov: int = 20, min_dep: int = 10,
                       sub_samp: int = 10, threads: int = 1) -> None:
    """
    Run complete microdiversity analysis.

    Args:
        output_dir: Base output directory
        ref_fasta: Path to reference FASTA file
        ref_genes: Path to reference genes file
        min_cov: Minimum coverage threshold
        min_dep: Minimum depth threshold
        sub_samp: Subsample size
        threads: Number of threads
    """
    base_path = os.path.normpath(output_dir)
    metapop_dir = os.path.join(base_path, "MetaPop")

    # Precompute Tajima's D parameters
    tajima_params = precompute_tajima_params(sub_samp + 10)

    # Get contig lengths
    contig_lengths = get_fasta_lengths(ref_fasta)

    # Process genic SNPs
    genic_snps_file = os.path.join(metapop_dir, "07.Cleaned_SNPs", "genic_snps.tsv")
    if os.path.exists(genic_snps_file):
        microdiv_data = process_genic_snps(genic_snps_file, sub_samp)

        # Calculate codon position summary
        codon_pos = calculate_codon_position_summary(microdiv_data)

        # Calculate contig-level microdiversity
        contig_microdiv = calculate_contig_level_microdiversity(
            microdiv_data, contig_lengths, tajima_params, sub_samp
        )

        # Save outputs
        output_path = os.path.join(metapop_dir, "10.Microdiversity")
        os.makedirs(output_path, exist_ok=True)

        if len(microdiv_data) > 0:
            # Drop sample_prop column to maintain column order compatibility with FST
            output_data = microdiv_data.drop(columns=['sample_prop'], errors='ignore')
            output_data.to_csv(
                os.path.join(output_path, "global_raw_microdiversity_data_snp_loci_only.tsv"),
                sep='\t', index=False
            )

        if len(codon_pos) > 0:
            codon_pos.to_csv(
                os.path.join(output_path, "global_codon_position_summary.tsv"),
                sep='\t', index=False
            )

        if len(contig_microdiv) > 0:
            contig_microdiv.to_csv(
                os.path.join(output_path, "global_contig_microdiversity.tsv"),
                sep='\t', index=False
            )

            # Check for suspicious codon position distribution
            total_third = codon_pos['third_pos'].sum()
            total_first = codon_pos['first_pos'].sum()
            total_second = codon_pos['second_pos'].sum()

            if total_third < total_second or total_third < total_first:
                print("Warning: More SNPs observed in first or second codon positions than third. "
                      "Results should be treated cautiously.")

        print("Microdiversity analysis complete.")
    else:
        print(f"Genic SNPs file not found: {genic_snps_file}")
