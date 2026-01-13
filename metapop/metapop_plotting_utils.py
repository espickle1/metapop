"""
MetaPop Plotting Utilities Module

Shared plotting utilities for MetaPop visualizations.
Provides matplotlib-based equivalents for R's ggplot2, cowplot, gggenes, etc.
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union
import warnings

# Try to import plotting libraries
try:
    import matplotlib
    matplotlib.use('Agg')  # Non-interactive backend for server/Colab
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import FancyArrow, FancyBboxPatch, Rectangle
    from matplotlib.collections import PatchCollection
    import matplotlib.colors as mcolors
    MATPLOTLIB_AVAILABLE = True
except ImportError:
    MATPLOTLIB_AVAILABLE = False
    warnings.warn("matplotlib not available. Plotting functions will not work.")

try:
    import seaborn as sns
    SEABORN_AVAILABLE = True
except ImportError:
    SEABORN_AVAILABLE = False


# =============================================================================
# Color Palette Functions (RColorBrewer equivalents)
# =============================================================================

def get_brewer_palette(name: str, n: int = 8) -> List[str]:
    """
    Get a ColorBrewer-style palette.

    Args:
        name: Palette name ('Spectral', 'RdYlBu', 'Set1', 'Set2', 'Set3', etc.)
        n: Number of colors

    Returns:
        List of hex color strings
    """
    if not MATPLOTLIB_AVAILABLE:
        return ['#000000'] * n

    palettes = {
        'Spectral': plt.cm.Spectral,
        'RdYlBu': plt.cm.RdYlBu,
        'RdYlGn': plt.cm.RdYlGn,
        'Set1': plt.cm.Set1,
        'Set2': plt.cm.Set2,
        'Set3': plt.cm.Set3,
        'Paired': plt.cm.Paired,
        'Pastel1': plt.cm.Pastel1,
        'Pastel2': plt.cm.Pastel2,
        'Dark2': plt.cm.Dark2,
        'Blues': plt.cm.Blues,
        'Reds': plt.cm.Reds,
        'Greens': plt.cm.Greens,
        'terrain': plt.cm.terrain,
    }

    cmap = palettes.get(name, plt.cm.viridis)
    colors = [mcolors.rgb2hex(cmap(i / (n - 1))) for i in range(n)]
    return colors


def get_pnps_colormap():
    """
    Get colormap for pN/pS ratio visualization.

    Returns:
        Colormap with blue (purifying) -> white (neutral) -> red (positive)
    """
    if not MATPLOTLIB_AVAILABLE:
        return None

    colors = ['#0000FF', '#FFFFFF', '#FF0000']  # Blue, white, red
    return mcolors.LinearSegmentedColormap.from_list('pnps', colors)


# =============================================================================
# Gene Arrow Plotting (gggenes equivalent)
# =============================================================================

def draw_gene_arrow(ax, start: int, end: int, strand: int, y: float = 0,
                    height: float = 0.8, color: str = 'gray',
                    edgecolor: str = 'black', label: str = None,
                    alpha: float = 1.0) -> None:
    """
    Draw a gene as an arrow on the given axes.

    Args:
        ax: Matplotlib axes
        start: Gene start position
        end: Gene end position
        strand: 1 for forward, -1 for reverse
        y: Y position for the arrow
        height: Height of the arrow
        color: Fill color
        edgecolor: Edge color
        label: Optional gene label
        alpha: Transparency
    """
    if not MATPLOTLIB_AVAILABLE:
        return

    width = abs(end - start)
    head_length = min(width * 0.15, 200)
    head_width = height * 0.6

    if strand == 1:  # Forward strand
        arrow = FancyArrow(
            start, y, width, 0,
            width=height,
            head_width=head_width + height,
            head_length=head_length,
            fc=color, ec=edgecolor, alpha=alpha,
            length_includes_head=True
        )
    else:  # Reverse strand
        arrow = FancyArrow(
            end, y, -width, 0,
            width=height,
            head_width=head_width + height,
            head_length=head_length,
            fc=color, ec=edgecolor, alpha=alpha,
            length_includes_head=True
        )

    ax.add_patch(arrow)

    if label:
        mid_x = (start + end) / 2
        ax.text(mid_x, y, label, ha='center', va='center',
                fontsize=6, rotation=0)


def plot_gene_track(ax, genes: pd.DataFrame, y: float = 0,
                    color_by: str = None, cmap: str = 'RdYlBu',
                    show_labels: bool = False) -> None:
    """
    Plot a track of genes as arrows.

    Args:
        ax: Matplotlib axes
        genes: DataFrame with columns: start, end, strand, [gene_name], [color_value]
        y: Y position for the track
        color_by: Column name to color genes by
        cmap: Colormap name for color_by values
        show_labels: Whether to show gene labels
    """
    if not MATPLOTLIB_AVAILABLE:
        return

    if color_by and color_by in genes.columns:
        values = genes[color_by].values
        norm = plt.Normalize(vmin=np.nanmin(values), vmax=np.nanmax(values))
        colormap = plt.cm.get_cmap(cmap)

    for _, gene in genes.iterrows():
        color = 'lightgray'
        if color_by and color_by in genes.columns:
            val = gene[color_by]
            if not pd.isna(val):
                color = colormap(norm(val))

        label = gene.get('gene_name', None) if show_labels else None

        draw_gene_arrow(
            ax,
            start=gene['start'],
            end=gene['end'],
            strand=gene.get('strand', 1),
            y=y,
            color=color,
            label=label
        )


# =============================================================================
# Depth/Coverage Plots
# =============================================================================

def plot_depth_coverage(ax, positions: np.ndarray, depths: np.ndarray,
                        color: str = 'steelblue', fill: bool = True,
                        alpha: float = 0.5, label: str = None) -> None:
    """
    Plot depth/coverage as a step plot.

    Args:
        ax: Matplotlib axes
        positions: Array of positions
        depths: Array of depth values
        color: Line/fill color
        fill: Whether to fill under the line
        alpha: Fill transparency
        label: Legend label
    """
    if not MATPLOTLIB_AVAILABLE:
        return

    ax.step(positions, depths, where='mid', color=color, label=label)
    if fill:
        ax.fill_between(positions, depths, step='mid', alpha=alpha, color=color)


# =============================================================================
# Heatmap Functions
# =============================================================================

def plot_clustered_heatmap(data: pd.DataFrame, output_path: str = None,
                           cmap: str = 'Spectral_r', figsize: Tuple = (10, 8),
                           show_rownames: bool = False, show_colnames: bool = True,
                           title: str = None) -> Optional[plt.Figure]:
    """
    Create a clustered heatmap similar to pheatmap in R.

    Args:
        data: DataFrame with data to plot
        output_path: Path to save figure (if None, returns figure)
        cmap: Colormap name
        figsize: Figure size
        show_rownames: Whether to show row names
        show_colnames: Whether to show column names
        title: Plot title

    Returns:
        Figure object if output_path is None
    """
    if not MATPLOTLIB_AVAILABLE or not SEABORN_AVAILABLE:
        print("matplotlib and seaborn required for heatmaps")
        return None

    # Remove all-zero rows
    data_nozero = data.loc[(data != 0).any(axis=1)]

    if len(data_nozero) == 0:
        print("No non-zero data to plot")
        return None

    g = sns.clustermap(
        data_nozero,
        cmap=cmap,
        figsize=figsize,
        yticklabels=show_rownames,
        xticklabels=show_colnames,
        dendrogram_ratio=0.1
    )

    if title:
        g.fig.suptitle(title, y=1.02)

    if output_path:
        g.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close(g.fig)
        return None
    else:
        return g.fig


# =============================================================================
# FST Heatmap
# =============================================================================

def plot_fst_heatmap(fst_matrix: pd.DataFrame, output_path: str = None,
                     title: str = 'FST Heatmap', figsize: Tuple = (8, 6)) -> None:
    """
    Plot FST values as a heatmap.

    Args:
        fst_matrix: Square matrix of FST values
        output_path: Path to save figure
        title: Plot title
        figsize: Figure size
    """
    if not MATPLOTLIB_AVAILABLE:
        return

    fig, ax = plt.subplots(figsize=figsize)

    if SEABORN_AVAILABLE:
        sns.heatmap(fst_matrix, ax=ax, cmap='YlOrRd', annot=True,
                   fmt='.3f', square=True, cbar_kws={'label': 'FST'})
    else:
        im = ax.imshow(fst_matrix.values, cmap='YlOrRd')
        plt.colorbar(im, ax=ax, label='FST')
        ax.set_xticks(range(len(fst_matrix.columns)))
        ax.set_xticklabels(fst_matrix.columns, rotation=45, ha='right')
        ax.set_yticks(range(len(fst_matrix.index)))
        ax.set_yticklabels(fst_matrix.index)

    ax.set_title(title)

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        plt.close()


# =============================================================================
# Scatter Plots with Regression
# =============================================================================

def plot_scatter_with_stats(ax, x: np.ndarray, y: np.ndarray,
                            xlabel: str = 'X', ylabel: str = 'Y',
                            title: str = None, show_mean: bool = True,
                            show_median: bool = True, color: str = 'gray') -> None:
    """
    Create scatter plot with mean/median lines.

    Args:
        ax: Matplotlib axes
        x: X values
        y: Y values
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        show_mean: Show mean horizontal line
        show_median: Show median horizontal line
        color: Point color
    """
    if not MATPLOTLIB_AVAILABLE:
        return

    ax.scatter(x, y, c=color, edgecolor='black', s=50, alpha=0.7)

    if show_mean:
        mean_val = np.nanmean(y)
        ax.axhline(mean_val, color='red', linestyle='-', label=f'Mean: {mean_val:.3f}')

    if show_median:
        median_val = np.nanmedian(y)
        ax.axhline(median_val, color='red', linestyle='--', label=f'Median: {median_val:.3f}')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    ax.legend(fontsize=8)


# =============================================================================
# Donut/Pie Charts
# =============================================================================

def plot_donut_chart(ax, values: List[float], labels: List[str],
                     colors: List[str] = None, title: str = None) -> None:
    """
    Create a donut chart.

    Args:
        ax: Matplotlib axes
        values: List of values
        labels: List of labels
        colors: List of colors
        title: Chart title
    """
    if not MATPLOTLIB_AVAILABLE:
        return

    if colors is None:
        colors = get_brewer_palette('Set2', len(values))

    wedges, texts, autotexts = ax.pie(
        values, labels=labels, colors=colors,
        autopct='%1.1f%%', pctdistance=0.75,
        wedgeprops=dict(width=0.5, edgecolor='white')
    )

    if title:
        ax.set_title(title)


# =============================================================================
# Bar Charts
# =============================================================================

def plot_stacked_bar(ax, data: pd.DataFrame, colors: List[str] = None,
                     xlabel: str = 'Samples', ylabel: str = 'Count',
                     title: str = None, show_legend: bool = True) -> None:
    """
    Create a stacked bar chart.

    Args:
        ax: Matplotlib axes
        data: DataFrame with categories as columns
        colors: List of colors for each category
        xlabel: X-axis label
        ylabel: Y-axis label
        title: Plot title
        show_legend: Whether to show legend
    """
    if not MATPLOTLIB_AVAILABLE:
        return

    if colors is None:
        colors = get_brewer_palette('Set2', len(data.columns))

    data.plot(kind='bar', stacked=True, ax=ax, color=colors, edgecolor='black')

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)
    if show_legend:
        ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.tick_params(axis='x', rotation=45)


# =============================================================================
# Multi-Panel Figures
# =============================================================================

def create_figure_grid(nrows: int, ncols: int,
                       figsize: Tuple = None) -> Tuple[plt.Figure, np.ndarray]:
    """
    Create a grid of subplots.

    Args:
        nrows: Number of rows
        ncols: Number of columns
        figsize: Figure size (auto-calculated if None)

    Returns:
        Tuple of (figure, axes array)
    """
    if not MATPLOTLIB_AVAILABLE:
        return None, None

    if figsize is None:
        figsize = (4 * ncols, 4 * nrows)

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)
    return fig, axes


def save_figure(fig: plt.Figure, output_path: str, dpi: int = 150) -> None:
    """
    Save figure to file.

    Args:
        fig: Matplotlib figure
        output_path: Output file path
        dpi: Resolution
    """
    if not MATPLOTLIB_AVAILABLE or fig is None:
        return

    fig.savefig(output_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
