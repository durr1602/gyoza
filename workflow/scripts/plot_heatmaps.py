from snakemake.script import snakemake
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use("Agg")  # Non-GUI backend for generating plots without display
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"

"""
This is a helper script to make sure we retrieve the necessary info
from dataframes of raw annotated read counts
to plot the corresponding heatmaps
"""

AA_LIST = "*PGCQNTSEDKHRWYFMLIVA"
AA_SORT = dict(zip(list(AA_LIST), list(range(0, len(AA_LIST)))))


def get_heatmaps(f, outpath, wtaa, level, label, plot_formats):
    # Extract time point
    t = outpath.split("_")[-1].split(".svg")[0]

    df = pd.read_csv(f, comment="#")

    idx = ["alt_aa"]
    if level == "codons":
        wtseq = df.loc[df.WT == True, "nt_seq"].values[0]
        wt_codons = [wtseq[i : i + 3] for i in range(0, len(wtseq), 3)]
        fitness = "readcount"
        idx.append("alt_codons")
        ccmap = plt.get_cmap("viridis")
    elif level == "aa":
        fitness = f"fitness_{t}"
        ccmap = sns.color_palette(
            "blend:#009B9E,#42B7B9,#A7D3D4,#F1F1F1,#E4C1D9,#D691C1,#C75DAB",  # CARTOColors Tropic
            as_cmap=True,
        )
    else:
        raise ValueError(
            f"Error.. Unrecognized level to plot heatmaps (supported: 'codons' or 'aa')"
        )

    singles = df[(df[f"Nham_{level}"] <= 1) & (df["aa_pos"] != "not-applicable")].copy()
    wide = singles.pivot(index=idx, columns="aa_pos", values=f"{fitness}")
    wide.sort_index(key=lambda x: x.map(AA_SORT), inplace=True)

    # Prepare graph space
    sns.set_theme(
        rc={
            "figure.figsize": (0.2 * len(wide.columns), 5),
            "font.size": 8,
            "legend.title_fontsize": 8,
            "legend.fontsize": 8,
            "axes.labelsize": 8,
            "axes.titlesize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
        },
        style="whitegrid",
    )
    f, ax = plt.subplots()

    # Color and mask for missing values
    ccmap.set_bad(".5")
    mask = pd.isnull(wide)

    # Draw heatmap
    ax = sns.heatmap(
        wide,
        mask=mask,
        cmap=ccmap,
        center=0,
        xticklabels=2,
        cbar_kws={"label": fitness, "shrink": 0.7},
    )

    ax.set_title(label)
    ax.set(xlabel=None, ylabel=None)
    ax.tick_params(axis="x", length=3)
    plt.yticks(rotation=0)

    cax = ax.figure.axes[-1]
    cax.tick_params(length=3, pad=2)

    # Build WT coordinates from dataframe if necessary
    if level == "codons":
        wtcoord = []
        for i, (aa, codon) in enumerate(zip(wtaa, wt_codons)):
            try:
                row_index = list(wide.index).index((aa, codon))
                wtcoord.append((i + 0.5, row_index + 0.5))
            except ValueError:
                # The (aa, codon) pair was not found in the data (skip)
                continue

    elif level == "aa":
        wtcoord = [(i + 0.5, list(AA_LIST).index(v) + 0.5) for i, v in enumerate(wtaa)]

    for o in wtcoord:
        ax.plot(o[0], o[1], marker="o", color="k", markersize=4)  # displays WT sequence

    # Export
    plt.tight_layout()
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]

    return


get_heatmaps(
    snakemake.input[0],
    snakemake.output[0],
    snakemake.params.wtaa,
    snakemake.params.level,
    snakemake.wildcards[0],
    [x for x in snakemake.config["plots"]["format"] if x != "svg"],
)
