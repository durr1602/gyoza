"""Module to plot heatmaps from mutants."""

from snakemake.script import snakemake
import pandas as pd
import pickle
import seaborn as sns
import matplotlib

matplotlib.use("Agg")  # Non-GUI backend for generating plots without display
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"


def get_heatmaps(f, meta, outpath, label, plot_formats):
    r"""Plot heatmap from dataframe, export in specified formats.
    
    Parameters
    ----------
    f : str
        Path to properly formatted dataframe.
    meta : str
        Path to serialized metadata. Should correspond to a dictionary with keys:
    
        * ``idx`` (**str** or **list of str**, index column headers)
        * ``fitness`` (**str**, column containing values to plot)
        * ``wt_coordinates`` (**list** of wild-type coordinates to mark on heatmap)
        * ``color_map``
        * ``vmax`` and ``vmin`` to set axes limits
    
    outpath : str
        Path to save heatmap as SVG (should end with ``.svg``).
    label : str
        Heatmap title
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    """
    # Load metadata
    with open(meta, "rb") as m:
        metadata = pickle.load(m)
    idx = metadata.get("idx")
    fitness = metadata.get("fitness")
    wtcoord = metadata.get("wt_coordinates")
    cmap = metadata.get("color_map")
    vmax = metadata.get("vmax")
    vmin = metadata.get("vmin")

    # Load dataframe
    df = pd.read_csv(f, index_col=idx)

    # Prepare graph space
    sns.set_theme(
        rc={
            "figure.figsize": (max(3, 0.2 * len(df.columns)), 0.2 * len(df.index)),
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

    # Mask for missing values
    mask = pd.isnull(df)

    # Draw heatmap
    ax = sns.heatmap(
        df,
        mask=mask,
        cmap=cmap,
        vmin=vmin,
        center=0,
        vmax=vmax,
        xticklabels=2,
        cbar_kws={"label": fitness, "shrink": 0.7},
    )

    ax.set_title(label)
    ax.set(xlabel=None, ylabel=None)
    ax.tick_params(axis="x", length=3)
    plt.yticks(rotation=0)

    cax = ax.figure.axes[-1]
    cax.tick_params(length=3, pad=2)

    # Display WT sequence
    for o in wtcoord:
        ax.plot(o[0], o[1], marker="o", color="k", markersize=4)

    # Export
    plt.tight_layout()
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]

    return


get_heatmaps(
    snakemake.input.heatmap_df,
    snakemake.input.heatmap_meta,
    snakemake.output[0],
    snakemake.wildcards[0],
    snakemake.params.plot_formats,
)
