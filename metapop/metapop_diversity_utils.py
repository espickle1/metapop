"""
MetaPop Diversity Utilities Module

Python implementations of vegan package diversity functions.
Replaces R's vegan, compositions, and related ecology packages.
"""

import numpy as np
import pandas as pd
from scipy.stats import gmean
from scipy.optimize import brentq
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.manifold import MDS
from typing import Tuple, Optional, Union, Dict


# =============================================================================
# Alpha Diversity Functions (vegan equivalents)
# =============================================================================

def specnumber(community_matrix: np.ndarray, axis: int = 1) -> np.ndarray:
    """
    Calculate species richness (number of species with abundance > 0).
    Equivalent to vegan::specnumber().

    Args:
        community_matrix: Samples x Species abundance matrix
        axis: Axis along which to count (1 = count species per sample)

    Returns:
        Array of richness values per sample
    """
    return np.sum(community_matrix > 0, axis=axis)


def diversity_shannon(abundances: np.ndarray) -> float:
    """
    Calculate Shannon's diversity index (H').
    Equivalent to vegan::diversity(x, index="shannon").

    H' = -sum(p_i * ln(p_i))

    Args:
        abundances: 1D array of species abundances

    Returns:
        Shannon's H value
    """
    abundances = abundances[abundances > 0]
    if len(abundances) == 0:
        return 0.0
    proportions = abundances / abundances.sum()
    return -np.sum(proportions * np.log(proportions))


def diversity_simpson(abundances: np.ndarray) -> float:
    """
    Calculate Simpson's diversity index (1 - D).
    Equivalent to vegan::diversity(x, index="simpson").

    D = sum(n_i * (n_i - 1)) / (N * (N - 1))
    Simpson = 1 - D

    Args:
        abundances: 1D array of species abundances

    Returns:
        Simpson's diversity value (1 - D)
    """
    n = abundances.sum()
    if n <= 1:
        return 0.0
    return 1 - np.sum(abundances * (abundances - 1)) / (n * (n - 1))


def diversity_invsimpson(abundances: np.ndarray) -> float:
    """
    Calculate Inverse Simpson's index (1/D).
    Equivalent to vegan::diversity(x, index="invsimpson").

    Args:
        abundances: 1D array of species abundances

    Returns:
        Inverse Simpson value
    """
    n = abundances.sum()
    if n <= 1:
        return 1.0
    d = np.sum(abundances * (abundances - 1)) / (n * (n - 1))
    if d == 0:
        return np.inf
    return 1 / d


def diversity(community_matrix: np.ndarray, index: str = "shannon") -> np.ndarray:
    """
    Calculate diversity indices for community matrix.
    Equivalent to vegan::diversity().

    Args:
        community_matrix: Samples x Species abundance matrix
        index: "shannon", "simpson", or "invsimpson"

    Returns:
        Array of diversity values per sample
    """
    func_map = {
        "shannon": diversity_shannon,
        "simpson": diversity_simpson,
        "invsimpson": diversity_invsimpson
    }

    if index not in func_map:
        raise ValueError(f"Unknown index: {index}. Use 'shannon', 'simpson', or 'invsimpson'")

    func = func_map[index]
    return np.array([func(row) for row in community_matrix])


def fisher_alpha(abundances: np.ndarray) -> float:
    """
    Calculate Fisher's alpha diversity index.
    Equivalent to vegan::fisher.alpha().

    Uses iterative solution for: S = alpha * ln(1 + N/alpha)
    where S = species count, N = total individuals

    Args:
        abundances: 1D array of species abundances

    Returns:
        Fisher's alpha value
    """
    S = np.sum(abundances > 0)  # Species count
    N = abundances.sum()  # Total individuals

    if S == 0 or N == 0:
        return np.nan

    # Solve: alpha * ln(1 + N/alpha) = S
    def equation(alpha):
        if alpha <= 0:
            return -S
        return alpha * np.log(1 + N / alpha) - S

    try:
        # Find alpha using Brent's method
        alpha = brentq(equation, 0.001, 10000)
        return alpha
    except ValueError:
        return np.nan


def chao1(abundances: np.ndarray) -> float:
    """
    Calculate Chao1 richness estimator.
    Equivalent to vegan::estimateR()["S.chao1"].

    Chao1 = S_obs + f1^2 / (2 * f2)
    where f1 = singletons, f2 = doubletons

    Args:
        abundances: 1D array of species abundances (integers)

    Returns:
        Chao1 estimated richness
    """
    S_obs = np.sum(abundances > 0)
    f1 = np.sum(abundances == 1)  # Singletons
    f2 = np.sum(abundances == 2)  # Doubletons

    if f2 > 0:
        return S_obs + (f1 ** 2) / (2 * f2)
    else:
        # Bias-corrected form when f2 = 0
        return S_obs + f1 * (f1 - 1) / 2


def ace(abundances: np.ndarray, rare_threshold: int = 10) -> float:
    """
    Calculate ACE (Abundance-based Coverage Estimator).
    Equivalent to vegan::estimateR()["S.ACE"].

    Args:
        abundances: 1D array of species abundances (integers)
        rare_threshold: Maximum abundance to be considered "rare"

    Returns:
        ACE estimated richness
    """
    abundances = np.asarray(abundances, dtype=int)

    S_rare = np.sum((abundances > 0) & (abundances <= rare_threshold))
    S_abund = np.sum(abundances > rare_threshold)

    if S_rare == 0:
        return S_abund

    # Count of individuals in rare species
    n_rare = abundances[abundances <= rare_threshold].sum()

    if n_rare == 0:
        return S_abund

    # Number of singletons among rare
    f1 = np.sum(abundances == 1)

    # Sample coverage estimate
    C_ace = 1 - f1 / n_rare

    if C_ace == 0:
        return np.inf

    # Calculate coefficient of variation
    sum_fi_i_i_minus_1 = 0
    for i in range(1, rare_threshold + 1):
        fi = np.sum(abundances == i)
        sum_fi_i_i_minus_1 += fi * i * (i - 1)

    gamma_sq = max(0, (S_rare / C_ace) * sum_fi_i_i_minus_1 / (n_rare * (n_rare - 1)) - 1)

    return S_abund + S_rare / C_ace + (f1 / C_ace) * gamma_sq


def estimateR(community_matrix: np.ndarray) -> pd.DataFrame:
    """
    Estimate species richness with Chao1 and ACE.
    Equivalent to vegan::estimateR().

    Args:
        community_matrix: Samples x Species abundance matrix

    Returns:
        DataFrame with Observed, Chao1, ACE for each sample
    """
    results = []
    for row in community_matrix:
        row = np.asarray(row, dtype=int)
        results.append({
            'Observed': np.sum(row > 0),
            'Chao1': chao1(row),
            'ACE': ace(row)
        })
    return pd.DataFrame(results)


def pielou_j(abundances: np.ndarray) -> float:
    """
    Calculate Pielou's evenness index (J').
    J' = H' / ln(S)

    Args:
        abundances: 1D array of species abundances

    Returns:
        Pielou's J value
    """
    S = np.sum(abundances > 0)
    if S <= 1:
        return np.nan
    H = diversity_shannon(abundances)
    return H / np.log(S)


# =============================================================================
# Beta Diversity Functions (vegan equivalents)
# =============================================================================

def vegdist(community_matrix: np.ndarray, method: str = "bray") -> np.ndarray:
    """
    Calculate pairwise dissimilarity matrix.
    Equivalent to vegan::vegdist().

    Args:
        community_matrix: Samples x Species abundance matrix
        method: "bray" (Bray-Curtis), "jaccard", or "euclidean"

    Returns:
        Square dissimilarity matrix
    """
    if method == "bray":
        # Bray-Curtis dissimilarity
        dist = pdist(community_matrix, metric=_bray_curtis)
    elif method == "jaccard":
        # Jaccard dissimilarity (presence/absence)
        binary = (community_matrix > 0).astype(float)
        dist = pdist(binary, metric='jaccard')
    elif method == "euclidean":
        dist = pdist(community_matrix, metric='euclidean')
    else:
        raise ValueError(f"Unknown method: {method}")

    return squareform(dist)


def _bray_curtis(u: np.ndarray, v: np.ndarray) -> float:
    """Calculate Bray-Curtis dissimilarity between two samples."""
    return np.sum(np.abs(u - v)) / np.sum(u + v) if np.sum(u + v) > 0 else 0


def clr_transform(data: np.ndarray, pseudocount: float = 1e-10) -> np.ndarray:
    """
    Centered log-ratio transformation.
    Equivalent to compositions::clr().

    CLR(x) = log(x / geometric_mean(x))

    Args:
        data: Abundance matrix (samples x species)
        pseudocount: Small value added to zeros

    Returns:
        CLR-transformed matrix
    """
    # Add pseudocount to zeros
    data = np.where(data == 0, pseudocount, data)

    # Calculate geometric mean per row
    gm = gmean(data, axis=1, keepdims=True)

    # CLR transformation
    return np.log(data / gm)


# =============================================================================
# Ordination Functions
# =============================================================================

def perform_pca(distance_matrix: np.ndarray, scale: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform PCA on distance matrix.
    Equivalent to prcomp().

    Args:
        distance_matrix: Square distance matrix
        scale: Whether to scale the data

    Returns:
        Tuple of (coordinates, explained_variance_ratio)
    """
    pca = PCA()
    if scale:
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler()
        distance_matrix = scaler.fit_transform(distance_matrix)

    coords = pca.fit_transform(distance_matrix)
    return coords, pca.explained_variance_ratio_


def perform_pcoa(distance_matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Perform Principal Coordinates Analysis (PCoA).
    Equivalent to vegan::capscale() with ~-1.

    Args:
        distance_matrix: Square distance matrix

    Returns:
        Tuple of (coordinates, eigenvalues)
    """
    n = distance_matrix.shape[0]

    # Double-center the distance matrix
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ (distance_matrix ** 2) @ H

    # Eigen decomposition
    eigenvalues, eigenvectors = np.linalg.eigh(B)

    # Sort by eigenvalue (descending)
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]

    # Keep only positive eigenvalues
    pos_idx = eigenvalues > 0
    eigenvalues = eigenvalues[pos_idx]
    eigenvectors = eigenvectors[:, pos_idx]

    # Calculate coordinates
    coords = eigenvectors * np.sqrt(eigenvalues)

    return coords, eigenvalues


def perform_nmds(distance_matrix: np.ndarray, n_components: int = 2,
                 max_iter: int = 300, random_state: int = 42) -> Tuple[np.ndarray, float]:
    """
    Perform Non-metric Multidimensional Scaling (NMDS).
    Equivalent to vegan::metaMDS().

    Args:
        distance_matrix: Square distance matrix
        n_components: Number of dimensions (k)
        max_iter: Maximum iterations
        random_state: Random seed

    Returns:
        Tuple of (coordinates, stress)
    """
    mds = MDS(n_components=n_components, dissimilarity='precomputed',
              max_iter=max_iter, random_state=random_state, metric=False)
    coords = mds.fit_transform(distance_matrix)
    stress = mds.stress_

    # Normalize stress similar to vegan
    # vegan uses stress = sqrt(stress / sum(d^2))
    d_sq_sum = np.sum(distance_matrix ** 2) / 2  # Upper triangle only
    if d_sq_sum > 0:
        stress = np.sqrt(stress / d_sq_sum)

    return coords, stress


# =============================================================================
# Tajima's D Parameter Functions
# =============================================================================

def tajima_d_parameters(sample_size: int) -> Dict[str, float]:
    """
    Calculate Tajima's D parameterization constants.

    Args:
        sample_size: Number of sequences in sample

    Returns:
        Dictionary with harmonic, harmonic2, b1, b2, c1, c2, e1, e2
    """
    if sample_size < 2:
        return {k: np.nan for k in ['harmonic', 'harmonic2', 'b1', 'b2', 'c1', 'c2', 'e1', 'e2']}

    n = sample_size
    # Harmonic numbers
    harmonic = np.sum(1.0 / np.arange(1, n))
    harmonic2 = np.sum(1.0 / (np.arange(1, n) ** 2))

    b1 = (n + 1) / (3 * (n - 1))
    b2 = (2 * (n ** 2 + n + 3)) / (9 * n * (n - 1))
    c1 = b1 - (1 / harmonic)
    c2 = b2 - ((n + 2) / (harmonic * n)) + harmonic2 / (harmonic ** 2)
    e1 = c1 / harmonic
    e2 = c2 / (harmonic ** 2 + harmonic2)

    return {
        'harmonic': harmonic,
        'harmonic2': harmonic2,
        'b1': b1,
        'b2': b2,
        'c1': c1,
        'c2': c2,
        'e1': e1,
        'e2': e2
    }


def calculate_tajimas_d(pi: float, num_snps: int, params: Dict[str, float]) -> float:
    """
    Calculate Tajima's D statistic.

    Args:
        pi: Nucleotide diversity
        num_snps: Number of segregating sites
        params: Dictionary from tajima_d_parameters()

    Returns:
        Tajima's D value
    """
    if num_snps == 0 or np.isnan(params['e1']):
        return np.nan

    theta = num_snps / params['harmonic']
    d_num = pi - theta
    d_denom = np.sqrt(params['e1'] * num_snps + params['e2'] * num_snps * (num_snps - 1))

    if d_denom == 0:
        return np.nan

    return d_num / d_denom


def calculate_alpha_diversity(community_matrix: np.ndarray) -> pd.DataFrame:
    """
    Calculate all alpha diversity metrics for a community matrix.

    Args:
        community_matrix: Samples x Species abundance matrix (samples in rows)

    Returns:
        DataFrame with diversity metrics for each sample
    """
    community_matrix = np.asarray(community_matrix, dtype=float)
    int_matrix = np.asarray(community_matrix, dtype=int)

    results = []
    for i, (row, int_row) in enumerate(zip(community_matrix, int_matrix)):
        results.append({
            'Richness': specnumber(row.reshape(1, -1))[0],
            'Shannons_H': diversity_shannon(row),
            'Simpson': diversity_simpson(row),
            'InvSimpson': diversity_invsimpson(row),
            'Fisher': fisher_alpha(int_row),
            'Pielous_J': pielou_j(row),
            'Chao1': chao1(int_row),
            'ACE': ace(int_row)
        })

    return pd.DataFrame(results)


def calculate_beta_diversity(community_matrix: np.ndarray,
                            methods: list = None) -> Dict[str, np.ndarray]:
    """
    Calculate multiple beta diversity distance matrices.

    Args:
        community_matrix: Samples x Species abundance matrix
        methods: List of methods ('bray', 'jaccard', 'clr_euclidean')

    Returns:
        Dictionary mapping method names to distance matrices
    """
    if methods is None:
        methods = ['bray', 'jaccard', 'clr_euclidean']

    results = {}

    if 'bray' in methods:
        results['bray'] = vegdist(community_matrix, 'bray')

    if 'jaccard' in methods:
        results['jaccard'] = vegdist(community_matrix, 'jaccard')

    if 'clr_euclidean' in methods:
        clr_data = clr_transform(community_matrix)
        results['clr_euclidean'] = vegdist(clr_data, 'euclidean')

    return results
