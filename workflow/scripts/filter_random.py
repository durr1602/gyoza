from snakemake.script import snakemake
import pandas as pd


def filter_random(df, max_Nham_aa, filtered_out, discarded_out):
    """
    Short script to filter mutants based on
    specified number of amino acid changes.
    Important to output the "unexpected" dataframe.
    """
    unfiltered_df = pd.read_csv(df, dtype={"mutation_aa_pos": str})
    unfiltered_df[unfiltered_df.Nham_aa <= max_Nham_aa].to_csv(
        filtered_out, index=False
    )
    unfiltered_df[unfiltered_df.Nham_aa > max_Nham_aa].to_csv(
        discarded_out, index=False
    )
    return


filter_random(
    snakemake.input[0],
    snakemake.params.Nham_aa_max,
    snakemake.output.filtered,
    snakemake.output.discarded,
)
