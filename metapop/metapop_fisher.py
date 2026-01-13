"""
MetaPop Fisher Exact Test Module

Performs Fisher's exact test for SNP linkage analysis and calculates
phi coefficient for effect size. Replaces MetaPop_Fisher_Exact.R.
"""

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from typing import Tuple, Optional


def linkage_test(all_ref: int, second_is_ref: int, first_is_ref: int, all_snp: int) -> Tuple[Optional[float], Optional[float]]:
    """
    Perform Fisher's exact test on a 2x2 contingency table and calculate phi coefficient.

    Args:
        all_ref: Count where both positions are reference
        second_is_ref: Count where first is SNP, second is reference
        first_is_ref: Count where first is reference, second is SNP
        all_snp: Count where both positions are SNPs

    Returns:
        Tuple of (p-value, phi coefficient), or (None, None) if invalid input
    """
    contingency_table = np.array([[all_ref, second_is_ref],
                                   [first_is_ref, all_snp]])

    if np.any(contingency_table < 0):
        return (None, None)

    # Fisher's exact test
    _, p_value = fisher_exact(contingency_table)

    # Phi coefficient calculation
    # phi = (ad - bc) / sqrt((a+b)(c+d)(a+c)(b+d))
    top = (all_ref * all_snp) - (second_is_ref * first_is_ref)

    bot_a = np.sqrt(all_ref + first_is_ref)
    bot_b = np.sqrt(all_ref + second_is_ref)
    bot_c = np.sqrt(all_snp + first_is_ref)
    bot_d = np.sqrt(all_snp + second_is_ref)
    bot_div = bot_a * bot_b * bot_c * bot_d

    if bot_div == 0:
        phi = np.nan
    else:
        phi = top / bot_div

    return (p_value, phi)


def process_linked_snps(input_file: str, output_file: Optional[str] = None) -> pd.DataFrame:
    """
    Process linked SNP data file, perform Fisher's exact test on each row,
    and calculate phi coefficients.

    Args:
        input_file: Path to tab-separated file with columns:
                   ref_count, ref_second, ref_first, snp_count
        output_file: Path to write results (defaults to overwriting input_file)

    Returns:
        DataFrame with fisher_p and phi_coef columns added
    """
    if output_file is None:
        output_file = input_file

    # Read the linked SNP data
    linked_data = pd.read_csv(input_file, sep='\t')

    # Apply linkage test to each row
    results = linked_data.apply(
        lambda row: linkage_test(
            row['ref_count'],
            row['ref_second'],
            row['ref_first'],
            row['snp_count']
        ),
        axis=1
    )

    # Unpack results into separate columns
    linked_data['fisher_p'] = [r[0] for r in results]
    linked_data['phi_coef'] = [r[1] for r in results]

    # Remove rows with NA p-values
    linked_data = linked_data.dropna(subset=['fisher_p'])

    # Write results
    linked_data.to_csv(output_file, sep='\t', index=False)

    return linked_data


def main():
    """Command-line interface matching original R script behavior."""
    import sys

    if len(sys.argv) < 2:
        print("Usage: python metapop_fisher.py <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    process_linked_snps(input_file)


if __name__ == "__main__":
    main()
