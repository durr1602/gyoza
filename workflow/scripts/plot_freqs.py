"""Plotting module for allele frequencies."""

from snakemake.script import snakemake
import pandas as pd
import seaborn as sns
import matplotlib

matplotlib.use("Agg")  # Non-GUI backend for generating plots without display
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"


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


def plot_allele_freq(df, outpath, plot_formats):
    r"""Plot distributions of allele frequencies for each sample group.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of allele frequencies. Should contain columns:

        * ``Sample attributes`` (**str**, sample group identifier)
        * ``frequency`` (**float**)
        * ``Timepoint`` (**str**)
        * ``Replicate`` (**str**, replicates are shown as split violins)

    outpath : str
        Path to save violin plot as SVG (should end with ``.svg``).
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    """
    labels = df["Sample attributes"].unique()
    timepoints = sorted(df.Timepoint.unique())
    mean_exp_freq = (
        df.groupby("Sample attributes")[["Mean_exp_freq"]].first().mean(axis=None)
    )
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


def get_allele_freq_plot(df_files, outpath, rc_level, reported_samples, plot_formats):
    r"""Aggregate data and plot distributions of allele frequencies.

    Parameters
    ----------
    df_files : list of str
        List of paths to CSV-formatted dataframes.
    outpath : str
        Path to save violin plot as SVG (should end with ``.svg``).
    rc_level : {"nt_seq", "barcode"}
        Level to which read counts are attributed.
    reported_samples : list of str
        List of samples to include in plot.
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    """
    df = concatenate_df(df_files)
    freq_per_seq = (
        df[df.Sample_name.isin(reported_samples)]
        .groupby(
            ["Sample attributes", rc_level, "Timepoint", "Replicate", "Mean_exp_freq"]
        )[["frequency"]]
        .first()
        .reset_index()
    )
    plot_allele_freq(freq_per_seq, outpath, plot_formats)
    return


get_allele_freq_plot(
    snakemake.input.freq_df,
    snakemake.output.rc_var_plot,
    snakemake.params.readcount_level,
    snakemake.params.reported_samples,
    snakemake.params.plot_formats,
)
