import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"
from upsetplot import from_indicators
from upsetplot import UpSet


def plot_stacked_barplot(df, outpath, exp_rc_per_sample, plot_formats):
    samples = df["Sample_name"].to_list()
    width = 0.5
    color_dict = dict(
        zip(
            ["OK", "Trimming", "Merging", "Aggregating", "Unexpected"],
            sns.color_palette("Spectral_r", 6)[1:],
        )
    )

    f, ax = plt.subplots(figsize=(20, 5))
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


def plot_rc_per_seq(
    df1, df2, outpath, sample_group, exp_rc_per_var, mean_exp_freq, plot_formats
):
    """
    Expects a dataframe of raw read counts and
    equivalent converted into read frequencies
    (normalized with sample depth).
    """
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))

    sns.histplot(df1, element="step", bins=50, common_norm=False, log_scale=10, ax=ax1)
    ax1.axvline(x=exp_rc_per_var, linestyle="--", color=".8")
    ax1.set(xlabel="Raw read count")

    sns.histplot(df2, element="step", bins=50, log_scale=10, common_norm=False, ax=ax2)
    ax2.axvline(x=10**mean_exp_freq, linestyle="--", color=".8")
    ax2.set(xlabel="Frequency")

    plt.subplots_adjust(top=0.9)
    plt.suptitle(f'Samples attributes: {" | ".join(sample_group)}')
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


# Define global vars that we are going to use in several plots
CSCORES = [1, 2, 3]
CSCORE_COLORS = ["green", "orange", "red"]


def plot_upset_TR(df, conditions, outpath, sample_group, plot_formats):
    """
    Plots overlap across time points and replicates
    in the form of an upsetplot, counting the number of unique sequences.
    Conditions argument correspond to columns in the dataframe..
    ..should be boolean and indicate whether or not the sequence is in the
    combination of time point/replicate
    """
    fig = plt.figure(figsize=(6, 6))
    upset_obj = UpSet(
        from_indicators(conditions, data=df),
        # show_percentages=True,
        show_counts=True,
        min_subset_size="1%",
        sort_by="cardinality",
        element_size=None,
        intersection_plot_elements=0,  # height of intersection barplot in matrix elements
        totals_plot_elements=2,  # width of totals barplot in matrix elements
    )

    upset_obj.add_stacked_bars(
        by="confidence_score", colors=dict(zip(CSCORES, CSCORE_COLORS)), elements=3
    )

    upset_obj.add_catplot(
        value="mean_input",
        kind="violin",
        cut=0,
        density_norm="count",
        log_scale=10,
        linewidth=0.5,
        elements=3,  # height in number of matrix elements
    )

    d = upset_obj.plot(
        fig=fig
    )  # Assigns all plots to a dictionary containing axes subplots - same keys as gridspec returned by upset_obj.make_grid()
    ax0 = d[
        "extra0"  # Key corresponding to 1st stacked barplot - confidence score ('intersections' = intersection barplot)
    ]
    ax1 = d["extra1"]  # Key corresponding to 1st catplot - read count for input samples

    ax0.set_ylabel("# Variants")
    ax0.legend(title="Confidence score")

    ax1.set_ylabel("Mean\nT0 freq.")

    plt.subplots_adjust(top=0.95)
    plt.suptitle(f'Samples attributes: {" | ".join(sample_group)}')

    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def plot_timepoint_corr(df, outpath, plot_formats):
    """
    Plots pairwise comparisons of selection coefficients
    to look at correlation between time points.
    """
    g = sns.pairplot(
        df,
        hue="confidence_score",
        hue_order=CSCORES,
        palette=dict(zip(CSCORES, CSCORE_COLORS)),
        plot_kws={"s": 8, "alpha": 0.2},
        height=1.5,
        corner=True,
    )
    g.tight_layout()
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return
