from snakemake.script import snakemake
import pandas as pd


def inject_WT(f, wt, outpath):
    """
    Short script to inject WT sequence
    in dataframe of raw read counts
    for annotation.
    """
    df = pd.read_csv(f)
    df["WT_seq"] = wt
    df.to_csv(outpath, index=False)
    return


inject_WT(
    snakemake.input[0],
    snakemake.params.wt,
    snakemake.output[0],
)
