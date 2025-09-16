"""Plotting module for aggregated data.

"""

from snakemake.script import snakemake
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
    r"""Opens and concatenates multiple dataframes.
    
    Parameters
    ----------
    df_files : list of str
        List of paths to CSV-formatted dataframes.
    
    Returns
    -------
    pandas.DataFrame
        Single concatenated dataframe.
    """
    list_df = []
    for f in df_files:
        list_df.append(pd.read_csv(f))
    df = pd.concat(list_df, ignore_index=True)
    return df


def plot_allele_freq(df, outpath, mean_exp_freq, plot_formats):
    r"""Plot distributions of allele frequencies for each sample group.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of allele frequencies. Should contain columns:
        + "Sample attributes" (str, sample group identifier)
        + "frequency" (float)
        + "Timepoint"
        + "Replicate" (replicates are shown as split violins)
    outpath : str
        Path to save violin plot as SVG (should end with ".svg").
    mean_exp_freq : float
        Log10 of average expected allele frequency.
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
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
    avg_label_length = sum(len(label) for label in labels) / len(labels)
    if avg_label_length > 20:  # wrap long labels
        labels = [
            "\n".join([a[i : i + 20] for i in range(0, len(a), 20)]) for a in labels
        ]
    g.set_xticklabels(labels, rotation=min(90, 4.5 * avg_label_length), ha="right")
    g.tight_layout()
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def get_allele_freq_plot(df_files, outpath, plot_formats):
    r"""Aggregate data and plot distributions of allele frequencies.
    
    Parameters
    ----------
    df_files : list of str
        List of paths to CSV-formatted dataframes.
    outpath : str
        Path to save violin plot as SVG (should end with ".svg").
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    """
    df = concatenate_df(df_files)
    mean_exp_freq = (
        df.groupby("Sample attributes")[["Mean_exp_freq"]].first().mean(axis=None)
    )
    plot_allele_freq(df, outpath, mean_exp_freq, plot_formats)
    return


def plot_scoeff_violin(df, outpath, plot_formats):
    r"""Plot distributions of functional impact scores for each sample group.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of functional impact scores. Should contain columns:
        + "Sample attributes" (str, sample group identifier)
        + "s" (float, functional impact score)
        + "Compared timepoints" (str)
        + "Replicate" (replicates are shown as split violins)
    outpath : str
        Path to save violin plot as SVG (should end with ".svg").
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
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
    avg_label_length = sum(len(label) for label in labels) / len(labels)
    if avg_label_length > 20:  # wrap long labels
        labels = [
            "\n".join([a[i : i + 20] for i in range(0, len(a), 20)]) for a in labels
        ]
    g.set_xticklabels(labels, rotation=min(90, 4.5 * avg_label_length), ha="right")
    g.tight_layout()
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def plot_impact_over_time(df, outpath, plot_formats):
    r"""Plot functional impact over time for each sample group.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of functional impact scores. Should contain columns:
        + "Sample attributes" (str, sample group identifier)
        + "s" (float, functional impact score)
        + "Compared timepoints" (str)
        + "Replicate" (replicates are shown as different markers)
        + "Nham_aa" (int, number of amino acid changes as different colors)
    outpath : str
        Path to save plot as SVG (should end with ".svg").
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    """
    g = sns.relplot(
        data=df.sort_values(by="Compared timepoints"),
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
    r"""Plot Spearman correlation between replicates as heatmaps.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of functional impact scores. Should contain columns:
        + "Sample attributes" (str, sample group identifier)
        + "Compared timepoints" (str)
        + `replicates`
    replicates : list of str
        List of replicates that should feature as columns in `df`.
        Spearman correlation coefficients are obtained for each pairwise comparison
        between replicates, for each sample group.
    outpath : str
        Path to save heatmap as SVG (should end with ".svg").
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
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
    r"""Plot correlation between first two replicates, for each sample group.
    
    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of functional impact scores. Should contain columns:
        + "Sample attributes" (str, sample group identifier)
        + "Compared timepoints" (str)
        + `replicates`
    replicates : list of str
        List of replicates that should feature as columns in `df`.
    outpath : str
        Path to save scatter plot as SVG (should end with ".svg").
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
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
    r"""Aggregate data and plot functional impact scores.
    
    Parameters
    ----------
    df_files : list of str
        List of paths to CSV-formatted dataframes.
    scoeff_plot_outpath : str
        Path to save violin plot of distributions of scores as SVG (should end with
        ".svg").
    s_time_plot_outpath : str
        Path to save plot of functional impact over time as SVG (should end with
        ".svg").
    heatmaps_outpath : str
        Path to save heatmap of Spearman correlation coefficients for pairwise
        comparisons of replicates as SVG (should end with ".svg").
    replicate_plot_outpath : str
        Path to save scatter plot of correlation between first two replicates (or
        against same replicate if only one) as SVG (should end with ".svg")
    plot_formats : list of str
        Formats other than SVG in which all plots should be saved.
    """
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
    snakemake.params.plot_formats,
)

get_s_plots(
    snakemake.input.aa_df,
    snakemake.output.scoeff_violin_plot,
    snakemake.output.s_through_time_plot,
    snakemake.output.replicates_heatmap_plot,
    snakemake.output.replicates_plot,
    snakemake.params.plot_formats,
)
