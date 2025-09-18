"""Plotting module for aggregated read count statistics."""

from snakemake.script import snakemake
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib

matplotlib.use("Agg")  # Non-GUI backend for generating plots without display
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"


def plot_stacked_barplot(df, outpath, exp_rc_per_sample, plot_formats):
    r"""Stacked bar plot of read count statistics.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of read count statistics. Should contain columns:
        + "Sample_name" (str, sample identifier)
        + "OK" (int, total number of reads that passed all filters)
        + "Trimming" (int, number of reads discarded at trimming stage)
        + "Merging" (int, number of reads discarded at merging step)
        + "Aggregating" (int, number of singletons)
        + "Unexpected" (int, number of reads for unexpected sequences)
        + "Contain_Ns" (int, number of reads discarded for containing Ns)
    outpath : str
        Path to save bar plot as SVG (should end with ".svg").
    exp_rc_per_sample : float
        Expected read count per sample.
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    """
    samples = df["Sample_name"].to_list()
    width = 0.5
    custom_palette = [
        (0.00784313725490196, 0.6196078431372549, 0.45098039215686275),
        (0.8705882352941177, 0.5607843137254902, 0.0196078431372549),
        (0.00392156862745098, 0.45098039215686275, 0.6980392156862745),
        (0.8352941176470589, 0.3686274509803922, 0.0),
        (0.8, 0.47058823529411764, 0.7372549019607844),
        (0.792156862745098, 0.5686274509803921, 0.3803921568627451),
    ]
    color_dict = dict(
        zip(
            ["OK", "Trimming", "Merging", "Aggregating", "Unexpected", "Contain_Ns"],
            custom_palette,
        )
    )

    f, ax = plt.subplots(figsize=(max(1, len(samples)), 5))
    bottom = np.zeros(len(df))

    for l in color_dict.keys():
        p = ax.bar(
            samples,
            df[l].values,
            width,
            label=l,
            bottom=bottom,
            color=color_dict[l],
        )
        bottom += df[l].values

    ax.set_yscale("log", base=10)
    ax.set(ylim=(1e4, 1e7), ylabel="Read count")

    ax.axhline(y=exp_rc_per_sample, linestyle="--", color=".8")
    ax.annotate("Aim", (-4, 1.1 * exp_rc_per_sample), color=".5")

    ax.xaxis.set_ticks(samples)
    ax.set_xticklabels(samples, rotation=45, ha="right")
    ax.legend(framealpha=0.9)

    plt.tight_layout()
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def plot_unexp_plot(df, outpath, plot_formats):
    r"""Plots distributions of read counts for unexpected variants.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of read counts for unexpected variants. Should contain columns:
        + "Sample_name" (str, sample identifier)
        + "readcount" (int)
    outpath : str
        Path to save plot as SVG (should end with ".svg").
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    """
    if df.empty:
        f, ax = plt.subplots(figsize=(4, 4))
        ax.text(0.5, 0.5, "No unexpected sequences", ha="center", va="center")
        ax.set_axis_off()  # hide axes
    else:
        sns.kdeplot(
            data=df,
            x="readcount",
            hue="Sample_name",
            common_norm=False,
            log_scale=True,
            legend=False,
        )
        plt.xlabel("Read count of unexpected variants")
        plt.savefig(outpath, format="svg", dpi=300)
        [
            plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
            for x in plot_formats
        ]
    return


def get_pooled_stats(
    sample_stats,
    sample_unexpected,
    csv_outpath,
    barplot_outpath,
    unexpplot_outpath,
    reported_samples,
    exp_rc_per_sample,
    plot_formats,
):
    r"""Aggregate and plot read count statistics.

    Parameters
    ----------
    sample_stats : list of str
        List of paths to CSV-formatted dataframes of read count statistics.
    sample_unexpected : list of str
        List of paths to CSV-formatted dataframes of unexpected variants.
    csv_outpath : str
        Path to save df with aggregated read count statistics for all samples.
    barplot_outpath : str
        Path to save stacked bar plot as SVG (should end with ".svg").
    unexpplot_outpath : str
        Path to save plot with distributions of read counts for unexpected
        variants as SVG (should end with ".svg")
    reported_samples : list of str
        List of samples to include in plots.
    exp_rc_per_sample : float
        Expected read count per sample.
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    """
    list_df = []
    for f in sample_stats:
        list_df.append(
            pd.read_csv(f)[
                [
                    "Sample_name",
                    "Total_raw_reads",
                    "Trimming",
                    "Merging",
                    "Aggregating",
                    "Contain_Ns",
                ]
            ]
        )

    stats_df = pd.concat(list_df, ignore_index=True)

    list_df = []
    for f in sample_unexpected:
        sample_unexp = pd.read_csv(f)[["Sample_name", "readcount"]]
        list_df.append(sample_unexp)

    unexp_all_seqs = pd.concat(list_df, ignore_index=True)

    plot_unexp_plot(
        unexp_all_seqs[unexp_all_seqs["Sample_name"].isin(reported_samples)],
        unexpplot_outpath,
        plot_formats,
    )

    unexp_df = unexp_all_seqs.groupby("Sample_name")[["readcount"]].sum()

    stacked_df = pd.merge(
        left=stats_df,
        right=unexp_df.rename(columns={"readcount": "Unexpected"}),
        on="Sample_name",
    )

    stacked_df["OK"] = stacked_df["Total_raw_reads"] - stacked_df[
        ["Trimming", "Merging", "Aggregating", "Unexpected", "Contain_Ns"]
    ].sum(axis=1)

    stacked_df.sort_index(inplace=True)

    stacked_df.to_csv(csv_outpath)

    plot_stacked_barplot(
        stacked_df[stacked_df["Sample_name"].isin(reported_samples)],
        barplot_outpath,
        exp_rc_per_sample,
        plot_formats,
    )

    return


get_pooled_stats(
    snakemake.input.read_stats,
    snakemake.input.unexpected,
    snakemake.output.all_stats,
    snakemake.output.rc_filter_plot,
    snakemake.output.unexp_rc_plot,
    snakemake.params.reported_samples,
    snakemake.params.exp_rc_per_sample,
    snakemake.params.plot_formats,
)
