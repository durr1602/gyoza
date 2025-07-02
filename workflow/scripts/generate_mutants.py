from snakemake.script import snakemake
from scripts.my_functions import load_codon_dic, get_single_double, get_nt_seq

# from pathlib import Path
import pandas as pd
import itertools


def generate_mutants(wtseq_path, outpath, codon_table, codon_mode):
    """
    Generates all single mutants and optionally all double mutants,
    based on the provided wild-type sequence
    """
    # Create directory if does not already exist
    # outdir.mkdir(parents=True, exist_ok=True)

    # Load codon dictionary
    codon_dic = load_codon_dic(codon_table)

    # Load dataframe with all wild-type sequences
    wt_df = pd.read_csv(wtseq_path)
    wt_df["nt_seq"] = wt_df["WT_seq"].str.upper()

    for mutated_seq in wt_df["Mutated_seq"].unique():
        sub_df = wt_df[wt_df["Mutated_seq"] == mutated_seq]

        if len(sub_df["WT_seq"].unique()) > 1:
            raise ValueError(
                f"Error.. For any mutated sequence, there should be a single possible WT sequence."
            )

        # Get all single mutants (and double mutants if necessary)
        mutants_df = get_single_double(sub_df, codon_dic, codon_mode)

        # Reconstruct nucleotide sequences from mutations
        mutants_df["nt_seq"] = list(
            itertools.starmap(
                get_nt_seq, zip(mutants_df["WT_seq"], mutants_df["mutations"])
            )
        )
        # Concatenate with WT
        # The mutations are dropped so we can annotate anything in downstream steps
        expmut_df = pd.concat(
            [wt_df, mutants_df.drop(["mutations"], axis=1)], ignore_index=True
        )

        # Convert to slightly lighter data types
        for x in expmut_df.columns:
            expmut_df[x] = expmut_df[x].astype("string")

        # Export dataframe
        expmut_df.to_csv(outpath, index=False, compression="gzip", chunksize=100_000)

    return


generate_mutants(
    snakemake.input[0],
    snakemake.output[0],
    snakemake.config["codon"]["table"],
    snakemake.config["codon"]["mode"],
)
