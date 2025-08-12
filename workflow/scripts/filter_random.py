from snakemake.script import snakemake
import pandas as pd


def filter_random(df, pre_discarded, max_Nham_aa, filtered_out, discarded_out):
    """
    Short script to filter mutants based on
    specified number of amino acid changes.
    The "unexpected" dataframe combines filtered out mutants
    **and** previously discarded indel mutants.
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
