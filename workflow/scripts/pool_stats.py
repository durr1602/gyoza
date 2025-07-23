from snakemake.script import snakemake

# from scripts.plotting_functions import plot_stacked_barplot, plot_unexp_plot
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib

matplotlib.use("Agg")  # Non-GUI backend for generating plots without display
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"


def plot_stacked_barplot(df, outpath, exp_rc_per_sample, plot_formats):
    samples = df["Sample_name"].to_list()
    width = 0.5
    color_dict = dict(
        zip(
            ["OK", "Trimming", "Merging", "Aggregating", "Unexpected"],
            sns.color_palette("Spectral_r", 6)[1:],
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
    # Empty plot if empty dataframe
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
    """
    Pool read count statistics from all samples.
    Merges info with unexpected sequences (seen in sequencing dataset but not expected or filtered out)
    """
    list_df = []
    for f in sample_stats:
        list_df.append(
            pd.read_csv(f)[
                ["Sample_name", "Total_raw_reads", "Trimming", "Merging", "Aggregating"]
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
        ["Trimming", "Merging", "Aggregating", "Unexpected"]
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
    float(snakemake.config["rc_aims"]["exp_rc_per_sample"]),
    [x for x in snakemake.config["plots"]["format"] if x != "svg"],
)
