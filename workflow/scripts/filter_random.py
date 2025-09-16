"""Module to filter mutants based on the number of amino acid changes.

"""

from snakemake.script import snakemake
import pandas as pd


def filter_random(df, pre_discarded, max_Nham_aa, filtered_out, discarded_out):
    r"""Discard mutants above specified number of amino acid changes.

    Parameters
    ----------
    df : str
        Path to CSV-formatted dataframe of mutants to filter.
        Should contain column "mutation_aa_pos" (casted as str) and "Nham_aa".
    pre_discarded : str
        Path to CSV-formatted dataframe of previously discarded indel mutants.
    max_Nham_aa : int
        Acceptable maximum number of amino acid changes.
    filtered_out : str
        Path to save output dataframe (mutants with number of amino acid changes
        <= `max_Nham_aa`).
    discarded_out : str
        Path to save output concatenated dataframe (previously discarded mutants
        and mutants with number of amino acid changes > `max_Nham_aa`).
    """
    unfiltered_df = pd.read_csv(df, dtype={"mutation_aa_pos": str})
    unfiltered_df[unfiltered_df.Nham_aa <= max_Nham_aa].to_csv(
        filtered_out, index=False
    )
    pre_discarded_df = pd.read_csv(pre_discarded)
    too_many_mut = unfiltered_df[unfiltered_df.Nham_aa > max_Nham_aa].copy()
    discarded = pd.concat([pre_discarded_df, too_many_mut], ignore_index=True)
    discarded.to_csv(discarded_out, index=False)
    return


filter_random(
    snakemake.input.no_indels,
    snakemake.input.indels,
    snakemake.params.Nham_aa_max,
    snakemake.output.filtered,
    snakemake.output.discarded,
)
