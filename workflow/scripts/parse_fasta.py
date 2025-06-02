from snakemake.script import snakemake
import pandas as pd


def get_read_counts(fasta_file, outpath, sample_name, wt, readcount_level):
    with open(fasta_file, "r") as file:
        entries = file.read().split(">")[1:]

    readcount = [int(x.split("size=")[1].split("\n")[0]) for x in entries]
    seqs = [x.split("size=")[1].split("\n", 1)[1].replace("\n", "") for x in entries]
    fasta_df = pd.DataFrame(
        list(zip(seqs, readcount)), columns=[readcount_level, "readcount"]
    )
    fasta_df["Sample_name"] = sample_name
    fasta_df["WT_seq"] = wt
    fasta_df.to_csv(outpath, index=False)

    return


get_read_counts(
    snakemake.input[0],
    snakemake.output[0],
    snakemake.wildcards.sample,
    snakemake.params[0],
    snakemake.config["barcode"]["rc_level"],
)
