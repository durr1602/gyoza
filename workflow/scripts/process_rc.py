from snakemake.script import snakemake

"""
from scripts.my_functions import get_confidence_score, get_mutation_type
from scripts.plotting_functions import (
    plot_rc_per_seq,
    plot_upset_TR,
    plot_timepoint_corr,
)
"""
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

CSCORES = [1, 2, 3]
CSCORE_COLORS = ["green", "orange", "red"]


def get_confidence_score(g, threshold):
    """
    Labels variants with a confidence score based on read count at T0
    """
    if (g >= threshold).all():  # Above threshold in all replicates
        return 1  # best confidence score
    elif (g >= threshold).any():  # Above threshold in at least 1 replicate
        return 2  # medium confidence score
    else:
        return 3  # low confidence score


def plot_rc_per_seq(
    df1, df2, outpath, sample_group, thresh, thresh_freq, plot_formats
):
    """
    Expects a dataframe of raw read counts and
    equivalent converted into read frequencies
    (normalized with sample depth).
    """
    fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2, figsize=(10, 4))

    sns.histplot(df1, element="step", bins=50, common_norm=False, log_scale=10, ax=ax1)
    ax1.axvline(x=thresh, linestyle="--", color=".8")
    ax1.set(xlabel="Raw read count")

    sns.histplot(df2, element="step", bins=50, log_scale=10, common_norm=False, ax=ax2)
    ax2.axvline(x=10**thresh_freq, linestyle="--", color=".8")
    ax2.set(xlabel="Frequency")

    plt.subplots_adjust(top=0.9)
    plt.suptitle(f"Samples attributes: {sample_group}")
    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


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
    plt.suptitle(f"Samples attributes: {sample_group}")

    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def plot_timepoint_corr(df, outpath, sample_group, plot_formats):
    """
    Plots pairwise comparisons of selection coefficients
    to look at correlation between time points.
    """
    # Check number of columns
    if len([x for x in df.columns if x != "confidence_score"]) <= 1:
        f, ax = plt.subplots(figsize=(max(4, 0.1 * len(sample_group)), 4))
        ax.text(0.5, 0.5, "Not enough time points to plot", ha="center", va="center")
        ax.set_axis_off()  # hide axes
    else:
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
        plt.subplots_adjust(top=0.9)
    plt.suptitle(f"Samples attributes: {sample_group}")

    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def get_selcoeffs(
    readcount_files,
    nbgen_path,
    outpath,
    avg_outpath,
    histplot_outpath,
    upsetplot_outpath,
    timepointsplot_outpath,
    freq_outpath,
    aa_df_outpath,
    sample_group,
    layout_path,
    sample_attributes,
    rc_level,
    barcode_attributes,
    rc_threshold,
    plot_formats,
):
    """
    Processes read counts and transforms them into functional impact scores.
    Samples sharing the same set of conditions (based on sample attributes)
    are pooled together. The resulting groups of samples are processed separately.
    """
    mutation_attributes = [
        "mutated_codon",
        "mutation_aa_pos",
        "mutation_alt_codons",
        "mutation_alt_aa",
        "mutation_type",
    ]

    sequence_attributes = [
        "nt_seq",
        "aa_seq",
        "Nham_nt",
        "Nham_aa",
        "Nham_codons",
        "aa_pos",
        "alt_aa",
    ]

    prot_seq_attributes = [
        "Nham_aa",
        "aa_seq",
        "aa_pos",
        "alt_aa",
    ]

    # Get back tuple from str wildcard
    sample_group_tuple = tuple(sample_group.split("__"))

    layout = pd.read_csv(
        layout_path,
        dtype={
            "Replicate": str,
        },
    ).set_index("Sample_name")

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
            },
        )
        groupdf = groupdf.merge(
            layout[["Timepoint", "Replicate"]],
            left_on="Sample_name",
            right_index=True,
        )
        df_list.append(groupdf.explode(["aa_pos", "alt_aa"]))

    df = pd.concat(df_list, ignore_index=True)

    # Add rows corresponding to variants not present in all replicates/time points
    df["TR"] = df["Timepoint"] + "_" + df["Replicate"]
    conditions = df.TR.unique()
    T0_conditions = [x for x in conditions if "T0" in x]

    upset = df.pivot_table(
        index=mutation_attributes + sequence_attributes + barcode_attributes,
        columns="TR",
        values="readcount",
        fill_value=0,
    ).reset_index(level=mutation_attributes + sequence_attributes + barcode_attributes)

    upset["confidence_score"] = upset[T0_conditions].apply(
        lambda row: get_confidence_score(row, rc_threshold), axis=1
    )
    sequence_attributes += ["confidence_score"]

    # Get total number of sequences
    tot_rc_level = upset[rc_level].nunique()

    # Determine how many sequences are labeled with lowest confidence_score
    low_conf_count = upset[upset["confidence_score"] == 3][rc_level].nunique()

    # Compute proportion
    low_conf_fraction = low_conf_count / tot_rc_level

    # Warn if more than 25%
    if low_conf_fraction > 0.25:
        cscore_statement = (
            f"Warning: More than 25% of your {rc_level}s are labeled with low confidence "
            f"(i.e., sequenced fewer than {rc_threshold} times in all replicates). "
            "Consider reviewing the config file and adjusting the rc_threshold parameter."
        )
        warnings.warn(cscore_statement)

    # Calculate frequencies
    freq = upset.copy()
    freq_conditions = [f"{x}_freq" for x in conditions]
    T0_freq = [x for x in freq_conditions if "T0" in x]

    if (freq[conditions].sum() == 0).any(axis=None):
        raise Exception(
            f"Oops.. at least one of your condition (and/or combination of time point and replicate) shows a null sample depth (no reads at all!)\n"
            f"Make sure your sample layout is OK. Unique combination of attributes should each have their T0 samples referencing the same FASTQ files."
        )

    freq[freq_conditions] = freq[conditions].add(1) / freq[conditions].sum()

    # Retrieve overall mean frequency corresponding to the specified read count threshold
    mean_thresh_freq = (
        np.log10(
            (rc_threshold + 1) / freq.groupby("nt_seq")[conditions].first().sum()
        )
    ).mean(axis=None)

    # Plot read count per sequence
    graph1df = freq.groupby("nt_seq")[conditions].first()
    graph2df = freq.groupby("nt_seq")[freq_conditions].first()
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
    upset_freq["mean_input"] = upset_freq[T0_freq].mean(axis=1)
    bool_conditions = [f"{x}_indicator" for x in conditions]
    upset_freq[bool_conditions] = upset_freq[conditions].astype(bool)
    upset_sub = (
        upset_freq.groupby("nt_seq")[
            bool_conditions + ["mean_input", "confidence_score"]
        ]
        .first()
        .rename(columns=dict(zip(bool_conditions, conditions)))
    )
    plot_upset_TR(upset_sub, conditions, upsetplot_outpath, sample_group, plot_formats)

    # Reshape dataframe
    longfreq = freq.melt(
        id_vars=mutation_attributes + sequence_attributes + barcode_attributes,
        value_vars=freq_conditions,
        var_name="TR_freq",
        value_name="frequency",
        ignore_index=False,
    ).reset_index(drop=True)
    longfreq["Timepoint"] = longfreq.TR_freq.apply(lambda x: x.split("_")[0])
    longfreq["Replicate"] = longfreq.TR_freq.apply(lambda x: x.split("_")[1])
    timepoints = sorted(longfreq.Timepoint.unique())

    # Output dataframe to plot distribution of allele frequencies
    longfreq_per_seq = (
        longfreq.groupby(
            sequence_attributes + barcode_attributes + ["Timepoint", "Replicate"]
        )[["frequency"]]
        .first()
        .reset_index()
    )
    # Save metadata in df for simplicity
    longfreq_per_seq["Sample attributes"] = sample_group
    longfreq_per_seq["Mean_exp_freq"] = mean_thresh_freq
    longfreq_per_seq.to_csv(freq_outpath, index=False)

    # Calculate Log2(fold-change) for every time point relative to T0
    # We cannot just go back to wide format since we need to melt Replicate --> pivot again
    freq_wide = longfreq.pivot(
        index=mutation_attributes
        + sequence_attributes
        + barcode_attributes
        + ["Replicate"],
        columns="Timepoint",
        values="frequency",
    )
    lfc_combinations = [(x, "T0") for x in timepoints[1:]]
    lfc_cols = [f'Lfc_{"_".join(x)}' for x in lfc_combinations]
    for i, v in enumerate(lfc_cols):
        freq_wide[v] = freq_wide.apply(
            lambda row: np.log2(
                row[lfc_combinations[i][0]] / row[lfc_combinations[i][1]]
            ),
            axis=1,
        )

    # Normalize with number of cellular generations
    nbgen_df = pd.read_csv(nbgen_path, dtype={"Replicate": str})
    # Select correct group
    nbgen_group = nbgen_df[
        nbgen_df[sample_attributes].apply(tuple, axis=1) == sample_group_tuple
    ]
    nbgen_wide = nbgen_group.pivot(
        index="Replicate", columns="Timepoint", values="Nb_gen"
    )
    nbgen_wide.columns = [f"{x}_gen" for x in nbgen_wide.columns]
    for i, x in enumerate(timepoints):
        if i in [0, 1]:
            pass
        else:
            nbgen_wide[f"cumul_{x}_gen"] = nbgen_wide[
                [f"{t}_gen" for t in timepoints[1:i] + [x]]
            ].sum(axis=1)
    for x in nbgen_wide.columns:
        if "cumul_" in x:
            nbgen_wide[x.split("cumul_")[1]] = nbgen_wide[x]
            nbgen_wide.drop(x, axis=1, inplace=True)

    lfc_wide = freq_wide.reset_index().merge(
        right=nbgen_wide.reset_index(), on="Replicate"
    )
    gen_cols = nbgen_wide.columns

    for x in list(zip(lfc_cols, gen_cols)):
        lfc_wide[x[0]] /= lfc_wide[x[1]]

    # Normalize with median of silent mutants
    syn = (
        lfc_wide[(lfc_wide.Nham_nt > 0) & (lfc_wide.Nham_aa == 0)]
        .groupby(["Replicate", "nt_seq"])[lfc_cols]
        .first()
        .reset_index()
    )
    mediansyn = syn.groupby("Replicate")[lfc_cols].median()
    mediansyn.columns = [x.replace("Lfc", "med") for x in mediansyn.columns]
    med_cols = mediansyn.columns

    # Calculate functional impact scores
    selcoeff_cols = [x.replace("Lfc", "s") for x in lfc_cols]
    s_wide = lfc_wide.merge(right=mediansyn.reset_index(), on="Replicate")
    for i, s in enumerate(selcoeff_cols):
        s_wide[s] = s_wide[lfc_cols[i]] - s_wide[med_cols[i]]

    # Retrieve wild-type protein sequence
    wtaa = s_wide.loc[s_wide.Nham_nt == 0, "aa_seq"].values[0]

    # Export full dataframe
    s_wide.to_csv(outpath)

    # Calculate median functional impact score (over synonymous codons), for each replicate separately
    median_df = (
        s_wide.groupby(["Replicate"] + prot_seq_attributes)[
            selcoeff_cols + ["confidence_score"]
        ]
        .agg(
            dict(
                zip(
                    selcoeff_cols + ["confidence_score"],
                    ["median"] * len(selcoeff_cols) + ["min"],
                )
            )
        )
        .reset_index(level=prot_seq_attributes)
    )

    # Plot correlation between time points
    dataset1_r1 = median_df.index[0]  # select first replicate only
    graphdf = median_df.loc[dataset1_r1].reset_index()[
        selcoeff_cols + ["confidence_score"]
    ]

    # Warn if many low confidence variants
    # Note: we do this here since there's the same number of variants for each replicate
    has_score = graphdf[selcoeff_cols].notna().any(axis=1)
    perc_counts = graphdf[has_score].groupby("confidence_score").size()
    perc_by_score = (
        perc_counts.div(perc_counts.sum()).to_frame("proportion").reset_index()
    )
    if (
        perc_by_score.loc[(perc_by_score.confidence_score == 1, "proportion")] < 0.75
    ).any():
        warnings.warn(
            f"Warning.. Group {sample_group} shows less than 75% of variants labeled with high confidence.\n"
            "Because only these variants are used to calculate a median score across replicates,"
            "you may want to double check the config file and adjust the rc_threshold parameter."
        )

    # Note: I decided to keep low confidence variants in this plot
    # but we filter right after
    plot_timepoint_corr(graphdf, timepointsplot_outpath, sample_group, plot_formats)

    # Filter only high confidence variants + reshape
    median_long = (
        median_df[median_df.confidence_score == 1]
        .melt(
            id_vars=prot_seq_attributes,
            value_vars=selcoeff_cols,
            var_name="Compared timepoints",
            value_name="s",
            ignore_index=False,
        )
        .reset_index()
    )
    # Rename column to keep only output time point (all are compared relative to T0)
    median_long["Compared timepoints"] = median_long["Compared timepoints"].apply(
        lambda x: x.split("_")[1]
    )
    # Save metadata in df for simplicity
    median_long["Sample attributes"] = sample_group

    # Output dataframe to plot more graphs (aggregating over sample groups)
    median_long.to_csv(aa_df_outpath, index=False)

    # Calculate median across replicates for high confidence variants
    avg_df = (
        median_df[median_df.confidence_score == 1]
        .groupby(prot_seq_attributes)[selcoeff_cols]
        .agg(
            [
                "median",
                lambda x: (
                    np.percentile(x.dropna(), 2.5) if len(x.dropna()) > 0 else np.nan
                ),
                lambda x: (
                    np.percentile(x.dropna(), 97.5) if len(x.dropna()) > 0 else np.nan
                ),
            ]
        )
    )

    # Rename columns
    cols_to_rename = [x for x in avg_df.columns if x not in prot_seq_attributes]
    new_names = []
    for c in cols_to_rename:
        if c[1] == "median":
            new_names.append(f"fitness_{c[0].split('_')[1]}")
        elif c[1] == "<lambda_0>":
            new_names.append(f"lower_err_{c[0].split('_')[1]}")
        elif c[1] == "<lambda_1>":
            new_names.append(f"upper_err_{c[0].split('_')[1]}")
        else:
            print("could not rename columns")

    avg_df.columns = new_names

    for x in new_names:
        if "lower_err_" in x:
            tp = x.split("lower_err_")[1]
            avg_df[x] = avg_df[f"fitness_{tp}"] - avg_df[x]
        elif "upper_err_" in x:
            tp = x.split("upper_err_")[1]
            avg_df[x] = avg_df[x] - avg_df[f"fitness_{tp}"]

    # Export dataframe with fitness and error values
    avg_df.reset_index().to_csv(avg_outpath, index=False)

    return


get_selcoeffs(
    snakemake.input.readcounts,
    snakemake.input.nbgen,
    snakemake.output.selcoeffs,
    snakemake.output.avg_scores,
    snakemake.output.hist_plot,
    snakemake.output.upset_plot,
    snakemake.output.timepoints_plot,
    snakemake.output.freq_df,
    snakemake.output.aa_df,
    snakemake.wildcards.group_key,
    snakemake.params.layout,
    snakemake.params.sample_attributes,
    snakemake.params.readcount_level,
    snakemake.params.barcode_attributes,
    snakemake.params.rc_threshold,
    snakemake.params.plot_formats,
)
