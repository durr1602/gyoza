"""Module to calculate and plot allele frequencies."""

from snakemake.script import snakemake
import pandas as pd
import json
import numpy as np
import seaborn as sns
import matplotlib

matplotlib.use("Agg")  # Non-GUI backend for generating plots without display
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"
from upsetplot import from_indicators
from upsetplot import UpSet
import warnings

MUTATION_ATTRIBUTES = [
    "mutated_codon",
    "mutation_aa_pos",
    "mutation_alt_codons",
    "mutation_alt_aa",
    "mutation_type",
]

SEQUENCE_ATTRIBUTES = [
    "nt_seq",
    "aa_seq",
    "Nham_nt",
    "Nham_aa",
    "Nham_codons",
    "aa_pos",
    "alt_aa",
    "wt_aa",
]

CSCORES = [1, 2, 3]
CSCORE_COLORS = ["green", "orange", "red"]


def load_read_counts(readcount_files, layout):
    """Load read counts for each sample, merge with layout and concatenate.

    Parameters
    ----------
    readcount_files : list of str
        List of paths to CSV-formatted dataframes of read counts per sample.
    layout : pandas.DataFrame
        Dataframe of sample layout.

    Returns
    -------
    pandas.DataFrame
    """
    df_list = []

    for f in readcount_files:
        groupdf = pd.read_csv(
            f,
            dtype={
                "WT": "boolean",  # Boolean type supports missing data
                "mutation_aa_pos": str,
            },
            converters={
                "aa_pos": json.loads,  # List
                "alt_aa": json.loads,  # List
                "wt_aa": json.loads,  # List
            },
        )
        groupdf = groupdf.merge(
            layout[["Timepoint", "Replicate"]],
            left_on="Sample_name",
            right_index=True,
        )
        df_list.append(groupdf.explode(["aa_pos", "alt_aa", "wt_aa"]))

    df = pd.concat(df_list, ignore_index=True)

    return df


def get_confidence_score(g, threshold):
    r"""Get confidence score based on read count at T0.

    Parameters
    ----------
    g : pandas.Series
        Read counts across replicates for a single sequence.
    threshold : int
        Read count threshold.

    Returns
    -------
    {1, 2, 3}
        Confidence score:

        * ``1``: high, read count above threshold in all replicates
        * ``2``: medium, read count above threshold in at least one replicate
        * ``3``: low, read count below threshold in all replicates
    """
    if (g >= threshold).all():
        return 1
    elif (g >= threshold).any():
        return 2
    else:
        return 3


def build_variant_matrix(df, sample_group, rc_level, barcode_attributes, rc_threshold):
    """Build matrix of variants across time points and replicates.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of read counts. Should contain columns:

        * ``Sample_name`` (**str**, unique sample identifier)
        * ``Timepoint`` (**str**)
        * ``Replicate`` (**str**)
        * ``WT`` (**boolean**, is WT or not)
        * ``mutation_aa_pos`` (**str**, position in the protein sequence
        at which the wild-type codon has been mutated
        * ``aa_pos`` (**list** of positions in the protein sequence
        at which the wild-type residue has been mutated
        * ``alt_aa`` (**list** of alternative amino acid residues at ``aa_pos``)
        * ``wt_aa`` (**list** of wild-type amino acid residues at ``aa_pos``)

    sample_group : str
        Sample group identifier.
        Should contain sample and screening attributes concatenated with ``__``.
    rc_level : {"nt_seq", "barcode"}
        Level to which read counts are attributed.
    barcode_attributes : list of str
        List of barcode attributes (includes `rc_level`).
    rc_threshold : int
        Threshold to label variants with a confidence score based on their
        read count at T0 across replicates.

    Returns
    -------
    df : pandas.DataFrame
        Matrix of variants
    TR_sample_dict : dict
        Dictionary mapping each combination of time points and replicates
        with the corresponding ``Sample_name`` value as in `df`.

    Warns
    -----
    UserWarning
        If less than 75% high confidence variants.

    """
    # Build mapping conditions <> samples
    df["TR"] = df["Timepoint"] + "_" + df["Replicate"]
    TR_sample_dict = (
        df[["TR", "Sample_name"]]
        .drop_duplicates()
        .set_index("TR")["Sample_name"]
        .to_dict()
    )
    T0_conditions = [x for x in TR_sample_dict.keys() if "T0" in x]

    # Add rows corresponding to variants not present in all replicates/time points
    upset = df.pivot_table(
        index=MUTATION_ATTRIBUTES + SEQUENCE_ATTRIBUTES + barcode_attributes,
        columns="TR",
        values="readcount",
        fill_value=0,
    ).reset_index(level=MUTATION_ATTRIBUTES + SEQUENCE_ATTRIBUTES + barcode_attributes)

    upset["confidence_score"] = upset[T0_conditions].apply(
        lambda row: get_confidence_score(row, rc_threshold), axis=1
    )

    # Get total number of sequences
    tot_rc_level = upset[rc_level].nunique()

    # Determine how many "high confidence" variants
    high_conf_count = upset[upset["confidence_score"] == 1][rc_level].nunique()

    # Compute proportion
    high_conf_fraction = high_conf_count / tot_rc_level

    # Warn if less than 75%
    if high_conf_fraction < 0.75:
        cscore_statement = (
            f"Warning: For group {sample_group}, less than 75% of your {rc_level}s are labeled with high confidence "
            f"(i.e., sequenced fewer than {rc_threshold} times in all replicates). "
            "Because only these variants are used to calculate a median score across replicates,"
            "consider reviewing the config file and adjusting the rc_threshold parameter."
        )
        warnings.warn(cscore_statement, UserWarning)

    return upset, TR_sample_dict


def plot_rc_per_seq(df1, df2, outpath, sample_group, thresh, thresh_freq, plot_formats):
    r"""Plot side-by-side distributions of read counts/frequencies.

    Parameters
    ----------
    df1 : pandas.DataFrame
        Dataframe of raw read counts.
    df2 : pandas.DataFrame
        Dataframe of allele frequencies (read counts normalized with sample depth).
    outpath : str
        Path to save plot as SVG (should end with ``.svg``).
    sample_group : str
        Sample group identifier
    thresh : int
        Read count threshold
    thresh_freq : float
        Log10 of `thresh` normalized with sample depth.
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    """
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))

    sns.histplot(df1, element="step", bins=50, common_norm=False, log_scale=10, ax=ax1)
    ax1.axvline(x=thresh, linestyle="--", color=".8")
    ax1.set(xlabel="Raw read count")

    sns.histplot(df2, element="step", bins=50, log_scale=10, common_norm=False, ax=ax2)
    ax2.axvline(x=10**thresh_freq, linestyle="--", color=".8")
    ax2.set(xlabel="Frequency")

    plt.subplots_adjust(top=0.9)
    plt.suptitle(f"{sample_group}")
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def plot_upset_TR(df, conditions, outpath, sample_group, plot_formats):
    r"""Plot overlap of unique sequences found across time points and replicates.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of sequences found across time points and replicates.
        Should contain columns:

        * `conditions` (**bool**, indicates if the sequence was in the combination of
          time point / replicate)
        * ``confidence_score``
        * ``mean_input`` (**float**, average read frequency at T0)

    conditions : list of str
        Columns in `df`, should refer to combinations of time point / replicate
    outpath : str
        Path to save upset plot as SVG (should end with ``.svg``).
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    """
    # Check number of conditions
    if len(conditions) < 2:
        f, ax = plt.subplots(figsize=(max(4, 0.1 * len(sample_group)), 4))
        ax.text(
            0.5, 0.5, "Not enough conditions\nto plot overlap", ha="center", va="center"
        )
        ax.set_axis_off()  # hide axes

    else:
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

        if df["mean_input"].values[0] != "not-applicable":
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

        ax0.set_ylabel("# Variants")
        ax0.legend(title="Confidence score")

        if "extra1" in d.keys():
            ax1 = d[
                "extra1"
            ]  # Key corresponding to 1st catplot - read count for input samples
            ax1.set_ylabel("Mean\nT0 freq.")

        plt.subplots_adjust(top=0.95)

    plt.suptitle(f"{sample_group}")

    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def get_frequencies(
    freq,
    TR_sample_dict,
    sample_group,
    barcode_attributes,
    rc_threshold,
    freq_outpath,
    histplot_outpath,
    upsetplot_outpath,
    reported_samples,
    plot_formats,
):
    """Calculate read frequencies for grouped samples.

    Parameters
    ----------
    freq : pandas.DataFrame
        Matrix of variants.
    TR_sample_dict : dict
        Dictionary mapping each combination of time points and replicates
        with the corresponding ``Sample_name`` value as in `df`.
    sample_group : str
        Sample group identifier.
        Should contain sample and screening attributes concatenated with ``__``.
    barcode_attributes : list of str
        List of barcode attributes (includes `rc_level`).
    rc_threshold : int
        Threshold to label variants with a confidence score based on their
        read count at T0 across replicates.
    freq_outpath : str
        Path to save output dataframe of allele frequencies for downstream
        processing (plots + calculation of functional impact scores).
    histplot_outpath : str
        Path to save plot with distributions of read counts/frequencies as SVG
        (should end with ``.svg``).
    upsetplot_outpath : str
        Path to save upset plot showing overlap of unique sequences found across
        time points and replicates, as SVG (should end with ``.svg``).
    reported_samples : list of str
        List of samples to include in plot.
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.

    """
    # Retrieve conditions
    conditions = list(TR_sample_dict)
    freq_conditions = [f"{x}_freq" for x in conditions]

    # Get conditions from samples marked for reporting
    reported_conditions = [
        k for k, v in TR_sample_dict.items() if v in reported_samples
    ]

    if (freq[conditions].sum() == 0).any(axis=None):
        raise Exception(
            f"Oops.. at least one of your condition (and/or combination of time point and replicate) shows a null sample depth (no reads at all!)\n"
            f"Make sure your sample layout is OK. Unique combination of attributes should each have their T0 samples referencing the same FASTQ files."
        )

    # Calculate frequencies
    freq[freq_conditions] = freq[conditions].add(1) / freq[conditions].sum()

    # Retrieve overall mean frequency corresponding to the specified read count threshold
    mean_thresh_freq = (
        np.log10((rc_threshold + 1) / freq.groupby("nt_seq")[conditions].first().sum())
    ).mean(axis=None)

    # Plot read count per sequence
    graph1df = freq.groupby("nt_seq")[reported_conditions].first()
    graph2df = freq.groupby("nt_seq")[
        [f"{x}_freq" for x in reported_conditions]
    ].first()
    plot_rc_per_seq(
        graph1df,
        graph2df,
        histplot_outpath,
        sample_group,
        rc_threshold,
        mean_thresh_freq,
        plot_formats,
    )

    # Plot overlap across time points and replicates
    upset_freq = freq.copy()
    T0_reported_conditions = [f"{x}_freq" for x in reported_conditions if "T0" in x]
    if T0_reported_conditions:
        upset_freq["mean_input"] = upset_freq[T0_reported_conditions].mean(axis=1)
    else:
        upset_freq["mean_input"] = "not-applicable"
    bool_conditions = [f"{x}_indicator" for x in reported_conditions]
    upset_freq[bool_conditions] = upset_freq[reported_conditions].astype(bool)
    upset_sub = (
        upset_freq.groupby("nt_seq")[
            bool_conditions + ["mean_input", "confidence_score"]
        ]
        .first()
        .rename(columns=dict(zip(bool_conditions, reported_conditions)))
    )
    plot_upset_TR(
        upset_sub, reported_conditions, upsetplot_outpath, sample_group, plot_formats
    )

    # Reshape dataframe
    longfreq = freq.melt(
        id_vars=MUTATION_ATTRIBUTES
        + SEQUENCE_ATTRIBUTES
        + barcode_attributes
        + ["confidence_score"],
        value_vars=freq_conditions,
        var_name="TR_freq",
        value_name="frequency",
        ignore_index=False,
    ).reset_index(drop=True)
    longfreq["Timepoint"] = longfreq.TR_freq.apply(lambda x: x.split("_")[0])
    longfreq["Replicate"] = longfreq.TR_freq.apply(lambda x: x.split("_")[1])
    longfreq["Mean_exp_freq"] = mean_thresh_freq
    longfreq["Sample attributes"] = sample_group
    longfreq["Sample_name"] = longfreq.TR_freq.apply(
        lambda x: TR_sample_dict.get(x.split("_freq")[0])
    )
    longfreq.to_csv(freq_outpath, index=False)

    return


def main(
    readcount_files,
    sample_group,
    layout,
    freq_outpath,
    histplot_outpath,
    upsetplot_outpath,
    reported_samples,
    plot_formats,
    rc_level,
    barcode_attributes,
    rc_threshold,
):
    """Convert read counts into frequencies and generate diagnostic plots.

    Parameters
    ----------
    readcount_files : list of str
        List of paths to CSV-formatted dataframes of read counts per sample.
    sample_group : str
        Sample group identifier.
        Should contain sample and screening attributes concatenated with ``__``.
    layout : pandas.DataFrame
        Dataframe of sample layout.
    freq_outpath : str
        Path to save output dataframe of allele frequencies for downstream
        processing (plots + calculation of functional impact scores).
    histplot_outpath : str
        Path to save plot with distributions of read counts/frequencies as SVG
        (should end with ``.svg``).
    upsetplot_outpath : str
        Path to save upset plot showing overlap of unique sequences found across
        time points and replicates, as SVG (should end with ``.svg``).
    reported_samples : list of str
        List of samples to include in plot.
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    rc_level : {"nt_seq", "barcode"}
        Level to which read counts are attributed.
    barcode_attributes : list of str
        List of barcode attributes (includes `rc_level`).
    rc_threshold : int
        Threshold to label variants with a confidence score based on their
        read count at T0 across replicates.

    """
    df = load_read_counts(readcount_files, layout)
    freq, TR_sample_dict = build_variant_matrix(
        df, sample_group, rc_level, barcode_attributes, rc_threshold
    )
    get_frequencies(
        freq,
        TR_sample_dict,
        sample_group,
        barcode_attributes,
        rc_threshold,
        freq_outpath,
        histplot_outpath,
        upsetplot_outpath,
        reported_samples,
        plot_formats,
    )


if __name__ == "__main__":
    main(
        snakemake.input.readcounts,
        snakemake.wildcards.group_key,
        snakemake.params.layout,
        snakemake.output.freq_df,
        snakemake.output.hist_plot,
        snakemake.output.upset_plot,
        snakemake.params.reported_samples,
        snakemake.params.plot_formats,
        snakemake.params.readcount_level,
        snakemake.params.barcode_attributes,
        snakemake.params.rc_threshold,
    )
