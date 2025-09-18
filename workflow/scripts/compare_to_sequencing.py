"""Module to filter **sequenced** expected mutants."""

from snakemake.script import snakemake
import pandas as pd


def get_observed_mutants(expmut, readcounts, outpath, unexp_outpath, readcount_level):
    r"""Compare expected and sequenced mutants, filter accordingly.

    Parameters
    ----------
    expmut : str
        Path to CSV-formatted dataframe of expected mutants.
        Should contain columns "Mutated_seq" (locus ID), and `readcount_level`.
    readcounts : str
        Path to CSV-formatted dataframe of sequenced mutants with read counts.
        Should contain column `readcount_level`.
    outpath : str
        Path to save output dataframe (mutants both expected and sequenced).
    unexp_outpath : str
        Path to save output dataframe (mutants sequenced but unexpected).
    readcount_level : {"nt_seq", "barcode"}
        Level to which read counts are attributed.
    """
    expmut_df = pd.read_csv(expmut)
    readcounts_df = pd.read_csv(readcounts)
    union_df = pd.merge(
        left=expmut_df,
        right=readcounts_df,
        how="outer",
        on=["Mutated_seq", readcount_level],
        indicator="Location",
    )

    overlap_df = union_df[union_df.Location == "both"].drop("Location", axis=1)
    if not overlap_df.empty:
        overlap_df.to_csv(outpath, index=False)
    else:
        empty_overlap = pd.DataFrame(columns=overlap_df.columns)
        empty_overlap.to_csv(outpath, index=False)

    # Sequences observed in the sequencing dataset, but not expected
    unexpected_df = union_df[union_df.Location == "right_only"].drop("Location", axis=1)
    if not unexpected_df.empty:
        unexpected_df.to_csv(unexp_outpath, index=False)
    else:
        empty_unexpected = pd.DataFrame(columns=readcounts_df.columns)
        empty_unexpected.to_csv(unexp_outpath, index=False)

    return


get_observed_mutants(
    snakemake.input.expected_mut_path,
    snakemake.input.readcount_path,
    snakemake.output.observed,
    snakemake.output.unexpected,
    snakemake.params.readcount_level,
)
