"""Module to format annotated read counts for downstream heatmap plotting."""

from snakemake.script import snakemake
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt


def get_heatmap_rc_data(f, outpath, meta_out, exp_rc):
    r"""Reshape dataframe of annotated read counts, extract and save metadata.

    Parameters
    ----------
    f : str
        Path to CSV-formatted dataframe of annotated read counts.
        Should contain columns "WT", "nt_seq", "aa_seq", "Nham_codons",
        "mutation_aa_pos", "mutation_alt_codons", "mutation_alt_aa" and "readcount"
    outpath : str
        Path to save reshaped dataframe.
    meta_out : str
        Path to save serialized metadata.
    exp_rc : float
        Expected read count per sample.
    """
    AA_LIST = "*PGCQNTSEDKHRWYFMLIVA"
    AA_SORT = dict(zip(list(AA_LIST), list(range(0, len(AA_LIST)))))

    df = pd.read_csv(f)

    # Retrieve wild-type
    wtseq = df.loc[df.WT == True, "nt_seq"].values[0]
    wtaa = df.loc[df.WT == True, "aa_seq"].values[0]
    wt_codons = [wtseq[i : i + 3] for i in range(0, len(wtseq), 3)]

    # Reshape dataframe
    filtered = (
        df[df.Nham_codons == 1]
        .groupby(["mutation_aa_pos", "mutation_alt_codons", "mutation_alt_aa"])[
            ["readcount"]
        ]
        .agg(Log10_readcount=("readcount", lambda x: np.log10(x.sum())))
        .reset_index()
    )
    filtered["mutation_aa_pos"] = filtered["mutation_aa_pos"].astype(int)
    wide = filtered.pivot(
        index=["mutation_alt_aa", "mutation_alt_codons"],
        columns="mutation_aa_pos",
        values="Log10_readcount",
    )
    wide.sort_index(key=lambda x: x.map(AA_SORT), inplace=True)

    # Export dataframe
    wide.to_csv(outpath)

    # Calculate WT coordinates
    wtcoord = []
    for i, (aa, codon) in enumerate(zip(wtaa, wt_codons)):
        try:
            row_index = list(wide.index).index((aa, codon))
            wtcoord.append((i + 0.5, row_index + 0.5))
        except ValueError:
            # The (aa, codon) pair was not found in the data (skip)
            continue

    # Set color map
    cmap = plt.get_cmap("viridis")
    cmap.set_bad(".5")

    # Assemble and save metadata
    meta = {
        "idx": ["mutation_alt_aa", "mutation_alt_codons"],
        "fitness": "Log10_readcount",
        "wt_coordinates": wtcoord,
        "color_map": cmap,
        "vmax": np.log10(exp_rc),
        "vmin": 0,
    }
    with open(meta_out, "wb") as m:
        pickle.dump(meta, m)

    return


get_heatmap_rc_data(
    snakemake.input[0],
    snakemake.output.heatmap_df,
    snakemake.output.heatmap_meta,
    snakemake.params.exp_rc,
)
