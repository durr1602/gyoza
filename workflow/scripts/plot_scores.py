from snakemake.script import snakemake

"""
from scripts.my_functions import concatenate_df
from scripts.plotting_functions import (
    plot_allele_freq,
    plot_scoeff_violin,
    plot_impact_over_time,
    plot_spearman_heatmaps,
    plot_replicate_scatter,
)
"""
import pandas as pd
from scipy import stats
import seaborn as sns
import matplotlib

matplotlib.use("Agg")  # Non-GUI backend for generating plots without display
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"

prot_seq_attributes = [
    "Nham_aa",
    "aa_seq",
    "aa_pos",
    "alt_aa",
]


def concatenate_df(df_files):
    """
    Takes list of paths to dataframes as input.
    Returns single concatenated dataframe.
    """
    list_df = []
    for f in df_files:
        list_df.append(pd.read_csv(f))
    df = pd.concat(list_df, ignore_index=True)
    return df


def plot_allele_freq(df, outpath, mean_exp_freq, plot_formats):
    """
    Plots distributions of allele frequencies for each sample group.
    Sample groups must be in column "Sample attributes".
    Replicates are shown as split violins.
    """
    labels = df["Sample attributes"].unique()
    timepoints = sorted(df.Timepoint.unique())
    g = sns.catplot(
        df,
        x="Sample attributes",
        y="frequency",
        row="Timepoint",
        row_order=timepoints,
        hue="Replicate",
        palette="hls",
        split=True,  # should work for more than 2 samples but might be ugly
        log_scale=10,
        kind="violin",
        cut=0,
        linewidth=1,
        inner="quart",
        height=2,
        aspect=0.8 * len(labels),
    )
    g.map(plt.axhline, y=10**mean_exp_freq, linestyle="--", color=".8")

    g.set_axis_labels("", "Frequency")
    g.set_titles(row_template="{row_name}")
    g.set_xticklabels(labels, rotation=45, ha="right")
    g.tight_layout()
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def get_allele_freq_plot(df_files, outpath, plot_formats):
    df = concatenate_df(df_files)
    mean_exp_freq = (
        df.groupby("Sample attributes")[["Mean_exp_freq"]].first().mean(axis=None)
    )
    plot_allele_freq(df, outpath, mean_exp_freq, plot_formats)
    return


def plot_scoeff_violin(df, outpath, plot_formats):
    """
    Plots distributions of functional impact scores for each sample group.
    Sample groups must be in column "Sample attributes".
    Replicates are shown as split violins.
    """
    labels = df["Sample attributes"].unique()
    timepoints = sorted(df["Compared timepoints"].unique())
    g = sns.catplot(
        df,
        x="Sample attributes",
        y="s",
        row="Compared timepoints",
        row_order=timepoints,
        hue="Replicate",
        palette="hls",
        split=True,  # should work for more than 2 samples but might be ugly
        kind="violin",
        cut=0,
        linewidth=1,
        inner="quart",
        height=2,
        aspect=0.8 * len(labels),
    )

    g.set_axis_labels("", "s")
    g.set_titles(row_template="{row_name}")
    g.set_xticklabels(labels, rotation=45, ha="right")
    g.tight_layout()
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def plot_impact_over_time(df, outpath, plot_formats):
    """
    Plots functional impact over time for each sample group.
    Sample groups must be in column "Sample attributes".
    Bands indicate spread (standard deviation) for each mutation type.
    """
    g = sns.relplot(
        data=df,
        x="Compared timepoints",
        y="s",
        col="Sample attributes",
        col_wrap=3,
        hue="Nham_aa",
        palette="hls",
        style="Replicate",
        kind="line",
        markers=True,
        errorbar="sd",
        height=1.5,
    )
    g.set(xlabel="")
    g.set_titles(col_template="{col_name}")
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def plot_spearman_heatmaps(df, replicates, outpath, plot_formats):
    """
    Plots heatmaps with Spearman correlation coefficient,
    for each pairwise comparison between replicates,
    for each sample group.
    Sample groups must be in column "Sample attributes".
    """
    labels = df["Sample attributes"].unique()
    timepoints = sorted(df["Compared timepoints"].unique())
    g = sns.FacetGrid(
        df,
        col="Sample attributes",
        col_order=labels,
        row="Compared timepoints",
        row_order=timepoints,
    )
    g.set_titles(col_template="{col_name}", row_template="{row_name}")

    for i, t in enumerate(timepoints):
        for j, l in enumerate(labels):
            yl = []
            for y in replicates:
                xl = []
                for x in replicates:
                    xval = df[
                        (df["Sample attributes"] == l)
                        & (df["Compared timepoints"] == t)
                    ][x].values
                    yval = df[
                        (df["Sample attributes"] == l)
                        & (df["Compared timepoints"] == t)
                    ][y].values
                    spearmanr, sp = stats.spearmanr(xval, yval)
                    xl.append(spearmanr)
                yl.append(xl)
            spearmanwide = pd.DataFrame(yl, columns=replicates, index=replicates)
            sns.heatmap(
                spearmanwide,
                vmin=0.5,
                vmax=1,
                cmap="viridis_r",
                annot=True,
                fmt=".2f",
                ax=g.axes[i][j],
            )
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def plot_replicate_scatter(df, replicates, outpath, plot_formats):
    """
    Plots correlation between first two replicates, for each sample group.
    replicates argument = list with values corresponding to each replicate.
    Sample groups must be in column "Sample attributes".
    """
    g = sns.lmplot(
        df,
        x=replicates[0],
        y=replicates[1],
        col="Sample attributes",
        col_wrap=3,
        hue="Compared timepoints",
        palette="mako",
        height=1.5,
        scatter_kws={"s": 8, "alpha": 0.2},
    )
    g.set_titles(col_template="{col_name}")
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def get_s_plots(
    df_files,
    scoeff_plot_outpath,
    s_time_plot_outpath,
    heatmaps_outpath,
    replicate_plot_outpath,
    plot_formats,
):
    df = concatenate_df(df_files)
    df["Replicate"] = df["Replicate"].astype(str)
    plot_scoeff_violin(df, scoeff_plot_outpath, plot_formats)
    plot_impact_over_time(df, s_time_plot_outpath, plot_formats)

    # Save list of replicates
    replicates = sorted(df.Replicate.unique())
    if len(replicates) == 1:
        firstTwoReplicates = replicates * 2
    else:
        firstTwoReplicates = replicates[:2]

    # Reshape dataframe
    repwide = df.pivot(
        index=prot_seq_attributes + ["Sample attributes", "Compared timepoints"],
        columns="Replicate",
        values="s",
    ).reset_index()

    plot_spearman_heatmaps(repwide, replicates, heatmaps_outpath, plot_formats)
    plot_replicate_scatter(
        repwide, firstTwoReplicates, replicate_plot_outpath, plot_formats
    )
    return


get_allele_freq_plot(
    snakemake.input.freq_df,
    snakemake.output.rc_var_plot,
    [x for x in snakemake.config["plots"]["format"] if x != "svg"],
)

get_s_plots(
    snakemake.input.aa_df,
    snakemake.output.scoeff_violin_plot,
    snakemake.output.s_through_time_plot,
    snakemake.output.replicates_heatmap_plot,
    snakemake.output.replicates_plot,
    [x for x in snakemake.config["plots"]["format"] if x != "svg"],
)
