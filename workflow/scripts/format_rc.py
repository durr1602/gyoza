"""Module to format annotated read counts for downstream heatmap plotting."""

from snakemake.script import snakemake
import pandas as pd
import numpy as np
import pickle
import matplotlib.pyplot as plt


def get_heatmap_rc_data(
    f, outpath, meta_out, wtseq, wtaa, pos_start, codon_table, exp_rc
):
    r"""Reshape dataframe of annotated read counts, extract and save metadata.

    Parameters
    ----------
    f : str
        Path to CSV-formatted dataframe of annotated read counts.
        Should contain columns:

        * ``nt_seq``
        * ``aa_seq``
        * ``Nham_codons``
        * ``mutation_aa_pos``
        * ``mutation_alt_codons``
        * ``mutation_alt_aa``
        * ``readcount``

    outpath : str
        Path to save reshaped dataframe.
    meta_out : str
        Path to save serialized metadata.
    wtseq : str
        Wild-type nucleotide sequence.
    wtaa : str
        Wild-type amino acid sequence.
    pos_start : int
        Starting position in the protein sequence.
    codon_table : pandas.DataFrame
        Codon table associating codons to amino acid residues.
        Should contain columns ``codon`` and ``aminoacid``, with
        all values in **upper case**.
    exp_rc : float
        Expected read count per sample.
    """
    AA_LIST = "*PGCQNTSEDKHRWYFMLIVA"
    AA_SORT = dict(zip(list(AA_LIST), list(range(0, len(AA_LIST)))))

    # First, build a sorted index from the codon table to get all codons
    gc_df = codon_table.assign(
        aa_rank=lambda d: d["aminoacid"].map(AA_SORT),
    ).sort_values(
        by=["aa_rank", "codon"],
        ascending=[True, True],
    )

    full_index = pd.MultiIndex.from_frame(
        gc_df[["aminoacid", "codon"]],
        names=["mutation_alt_aa", "mutation_alt_codons"],
    )

    # Import read counts
    df = pd.read_csv(f)

    # Retrieve wild-type coordinates + positions
    wt_codons = [wtseq[i : i + 3] for i in range(0, len(wtseq), 3)]
    positions = np.arange(pos_start, pos_start + len(wtaa))

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

    # Reindex to add all codons and all positions (already sorted)
    wide = wide.reindex(index=full_index, columns=positions)

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
    snakemake.input.readcounts,
    snakemake.output.heatmap_df,
    snakemake.output.heatmap_meta,
    snakemake.params.wt["nt"],
    snakemake.params.wt["aa"],
    snakemake.params.wt["pos_start"],
    snakemake.params.codon_table,
    snakemake.params.exp_rc,
)
