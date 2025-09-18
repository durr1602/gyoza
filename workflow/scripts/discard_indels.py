"""Module to discard indels."""

from snakemake.script import snakemake
import pandas as pd


def inject_WT(f, wt, outpath, indel_outpath):
    r"""Inject WT sequence in dataframe of raw read counts, discard indel mutants.

    Parameters
    ----------
    f : str
        Path to CSV-formatted dataframe of raw read counts.
        Should contain column "nt_seq".
    wt : str
        Wild-type DNA sequence (with bases either A, C, G or T).
    outpath : str
        Path to save output dataframe (non-indel mutants only),
        with column headers only if empty.
    indel_outpath : str
        Path to save output dataframe (indels only),
        with column headers only if empty.
    """
    df = pd.read_csv(f)
    df["WT_seq"] = wt

    # Identify indels
    is_indel = df["nt_seq"].str.len() != df["WT_seq"].str.len()
    df_indels = df[is_indel].copy()
    df_valid = df[~is_indel].copy()

    # Ensure dataframe with indels is created even if empty
    if df_indels.empty:
        df_indels = pd.DataFrame(columns=df.columns)
    else:
        df_indels.to_csv(indel_outpath, index=False)

    df_valid.to_csv(outpath, index=False)

    return


inject_WT(
    snakemake.input[0],
    snakemake.params.wt,
    snakemake.output.observed,
    snakemake.output.indels,
)
