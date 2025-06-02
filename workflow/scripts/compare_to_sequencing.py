from snakemake.script import snakemake
import pandas as pd


def get_observed_mutants(expmut, readcounts, outpath, unexp_outpath, readcount_level):
    """
    Compares the dataframe of mutants observed in the sequencing data (with respective read count)
    to the list of expected mutants (either generated upstream or user-provided).
    Outputs the overlap = 1 file per sample containing only the mutants both expected and observed
    """
    expmut_df = pd.read_csv(expmut)
    readcounts_df = pd.read_csv(readcounts)
    union_df = pd.merge(
        left=expmut_df,
        right=readcounts_df,
        how="outer",
        on=["WT_seq", readcount_level],
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
    str(snakemake.input["expected_mut_path"]),
    str(snakemake.input["readcount_path"]),
    snakemake.output.observed,
    snakemake.output.unexpected,
    snakemake.config["barcode"]["rc_level"],
)
