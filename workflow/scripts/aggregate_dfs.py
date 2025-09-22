"""Module to concatenate final dataframes of scores."""

from snakemake.script import snakemake
import pandas as pd


def concatenate_dfs(input_dfs, output_df):
    r"""Concatenate a list of dataframes into a single one.
    
    Parameters
    ----------
    input_dfs : list of str
        List of paths to CSV-formatted dataframes.
    output_df : str
        Path to save output dataframe.
    """
    df_list = []
    for f in input_dfs:
        df_list.append(pd.read_csv(f))
    df = pd.concat(df_list, ignore_index=True)
    df.to_csv(output_df, index=None)

    return


concatenate_dfs(
    snakemake.input.all_df,
    snakemake.output.selcoeffs,
)

concatenate_dfs(
    snakemake.input.avg_df,
    snakemake.output.avg_scores,
)
