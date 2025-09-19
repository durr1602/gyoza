"""Module to convert fasta file to dataframe of read counts."""

from snakemake.script import snakemake
import pandas as pd


def get_read_counts(
    fasta_file, outpath, log_file, mutated_seq, sample_name, readcount_level
):
    r"""Parse fasta file, convert to dataframe of read counts, filter and save.
    
    Parameters
    ----------
    fasta_file : str
        Path to fasta file. Sequence headers should contain ``size={readcount}``.
    outpath : str
        Path to save output dataframe of read counts (only sequences with no Ns).
    log_file : str
        Path to custom log file to save the number of sequences discarded because
        they contained Ns.
    mutated_seq : str
        Locus identifier
    sample_name : str
        Sample identifier
    readcount_level : {"nt_seq", "barcode"}
        Level to which read counts are attributed.
    """
    with open(fasta_file, "r") as file:
        entries = file.read().split(">")[1:]

    readcount = [int(x.split("size=")[1].split("\n")[0]) for x in entries]
    seqs = [x.split("size=")[1].split("\n", 1)[1].replace("\n", "") for x in entries]
    fasta_df = pd.DataFrame(
        list(zip(seqs, readcount)), columns=[readcount_level, "readcount"]
    )

    # Add metadata
    fasta_df["Mutated_seq"] = mutated_seq
    fasta_df["Sample_name"] = sample_name

    # Discard sequences with Ns and log number
    withNs = fasta_df[readcount_level].str.contains("N")
    filtered_df = fasta_df[~withNs].copy()
    n_discarded = len(fasta_df) - len(filtered_df)
    with open(log_file, "w") as file:
        file.write(f"{sample_name},{n_discarded}")

    # Export sequences with no Ns
    filtered_df.to_csv(outpath, index=False)

    return


get_read_counts(
    snakemake.input[0],
    snakemake.output.readcounts,
    snakemake.log[0],
    snakemake.params.mutseq,
    snakemake.wildcards.sample,
    snakemake.params.readcount_level,
)
