from snakemake.script import snakemake
from scripts.my_functions import get_confidence_score, get_mutation_type
from scripts.plotting_functions import (
    plot_rc_per_seq,
    plot_upset_TR,
    plot_timepoint_corr,
)
import pandas as pd
import numpy as np
import warnings


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
    barcode_attributes,
    rc_level,
    rc_threshold,
    exp_rc_per_var,
    plot_formats,
):
    """
    Processes read counts and transforms them into functional impact scores.
    Samples sharing the same set of conditions (based on sample attributes)
    are pooled together. The resulting groups of samples are processed separately.
    """
    mutation_attributes = [
        "mutated_codon",
        "pos",
        "alt_codons",
        "alt_aa",
        "aa_pos",
        "nt_seq",
        "aa_seq",
        "Nham_nt",
        "Nham_aa",
        "Nham_codons",
    ]

    mutation_attributes_aa = [
        "mutated_codon",
        "aa_pos",
        "alt_aa",
        "Nham_aa",
        "mutation_type",
    ]

    layout = pd.read_csv(layout_path).set_index("Sample_name")

    df_list = []

    for f in readcount_files:
        groupdf = pd.read_csv(
            f,
            dtype={
                "WT": "boolean",  # Boolean type supports missing data
                "pos": str,
            },
        )
        groupdf = groupdf.merge(
            layout[["Pos_start", "Timepoint", "Replicate"]],
            left_on="Sample_name",
            right_index=True,
        )
        df_list.append(groupdf)

    df = pd.concat(df_list, ignore_index=True)

    # Add position offset
    position_offset = df.Pos_start.unique().item()
    df["aa_pos"] = df.pos.apply(
        lambda x: (
            int(x) + position_offset if x != "not-applicable" else "not-applicable"
        )
    )

    # Add rows corresponding to variants not present in all replicates/time points
    df["TR"] = df["Timepoint"] + "_" + df["Replicate"]
    conditions = df.TR.unique()
    T0_conditions = [x for x in conditions if "T0" in x]

    upset = df.pivot_table(
        index=mutation_attributes + barcode_attributes,
        columns="TR",
        values="readcount",
        fill_value=0,
    ).reset_index(level=mutation_attributes + barcode_attributes)

    upset["confidence_score"] = upset[T0_conditions].apply(
        lambda row: get_confidence_score(row, rc_threshold), axis=1
    )
    mutation_attributes += ["confidence_score"]

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

    # Retrieve overall mean frequency corresponding to the expected read count per variant
    mean_exp_freq = (np.log10((exp_rc_per_var + 1) / freq[conditions].sum())).mean(
        axis=None
    )

    # Plot read count per sequence
    graph1df = freq[conditions]
    graph2df = freq[freq_conditions]
    plot_rc_per_seq(
        graph1df,
        graph2df,
        histplot_outpath,
        sample_group,
        exp_rc_per_var,
        mean_exp_freq,
        plot_formats,
    )

    # Plot overlap across time points and replicates
    freq["mean_input"] = freq[T0_freq].mean(axis=1)
    bool_conditions = [f"{x}_indicator" for x in conditions]
    freq[bool_conditions] = freq[conditions].astype(bool)
    upset_sub = freq.drop(conditions, axis=1).rename(
        columns=dict(zip(bool_conditions, conditions))
    )
    plot_upset_TR(upset_sub, conditions, upsetplot_outpath, sample_group, plot_formats)

    # Output dataframe to plot distribution of allele frequencies
    longfreq = freq.melt(
        id_vars=mutation_attributes + barcode_attributes,
        value_vars=freq_conditions,
        var_name="TR_freq",
        value_name="frequency",
        ignore_index=False,
    ).reset_index(drop=True)
    longfreq["Timepoint"] = longfreq.TR_freq.apply(lambda x: x.split("_")[0])
    longfreq["Replicate"] = longfreq.TR_freq.apply(lambda x: x.split("_")[1])
    timepoints = sorted(longfreq.Timepoint.unique())
    longfreq.to_csv(freq_outpath, index=False)

    # Annotate with mutation type
    longfreq["mutation_type"] = longfreq.apply(
        lambda row: get_mutation_type(row.Nham_aa, row.alt_aa), axis=1
    )
    mutation_attributes += ["mutation_type"]

    # Calculate Log2(fold-change) for every time point relative to T0
    freq_wide = longfreq.pivot(
        index=mutation_attributes + ["Replicate"] + barcode_attributes,
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
    nbgen_df = pd.read_csv(nbgen_path)
    nbgen_wide = nbgen_df.pivot(index="Replicate", columns="Timepoint", values="Nb_gen")
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
    syn = lfc_wide[(lfc_wide.Nham_nt > 0) & (lfc_wide.Nham_aa == 0)][
        ["Replicate"] + lfc_cols
    ]
    mediansyn = syn.groupby("Replicate")[lfc_cols].median()
    mediansyn.columns = [x.replace("Lfc", "med") for x in mediansyn.columns]
    med_cols = mediansyn.columns

    # Calculate functional impact scores
    selcoeff_cols = [x.replace("Lfc", "s") for x in lfc_cols]
    s_wide = lfc_wide.merge(right=mediansyn.reset_index(), on="Replicate")
    for i, s in enumerate(selcoeff_cols):
        s_wide[s] = s_wide[lfc_cols[i]] - s_wide[med_cols[i]]

    ### Repeat WT at every position
    # Select WT nucleotide sequence(s)
    WTdf = s_wide[s_wide.Nham_nt == 0]

    # Get length of protein sequence
    WTdf["len_aa"] = WTdf.aa_seq.apply(lambda x: len(x))

    # Create list of positions for every sequence
    WTdf["pos"] = WTdf.len_aa.apply(lambda x: np.arange(x))

    # Same with WT codons (list of codons at every matching position)
    WTdf["alt_codons"] = WTdf.nt_seq.apply(
        lambda x: [x[i : i + 3] for i in range(0, len(x), 3)]
    )

    # Same with WT amino acid
    WTdf["alt_aa"] = WTdf.aa_seq.apply(lambda x: [y for y in x])

    # Then we use explode to turn horizontal lists into rows with matching values for all 3 columns
    WTdf = WTdf.explode(["pos", "alt_codons", "alt_aa"])

    # And finally add the position offset (position in the full protein sequence)
    WTdf["aa_pos"] = WTdf["pos"] + position_offset

    # Get non-WT
    # In this step we need to cast the dtype of pos and aa_pos
    # which we could not do before because the WT rows feature string values ("non-applicable")
    nonWT = s_wide[s_wide.Nham_nt > 0]
    nonWT[["pos", "aa_pos"]] = nonWT[["pos", "aa_pos"]].astype(int)

    # Concatenate to get full dataframe of unaggregated scores
    allpos_df = pd.concat([WTdf, nonWT], ignore_index=True)

    # Export
    allpos_df.to_csv(outpath)

    # Calculate median functional impact score (over synonymous mutants)
    median_df = (
        allpos_df.groupby(["Replicate", "mutated_codon", "aa_pos", "alt_aa"])[
            selcoeff_cols + ["confidence_score", "Nham_aa", "mutation_type"]
        ]
        .agg(
            dict(
                zip(
                    selcoeff_cols + ["confidence_score", "Nham_aa", "mutation_type"],
                    ["median"] * len(selcoeff_cols) + ["min", "first", "first"],
                )
            )
        )
        .reset_index(level=["mutated_codon", "aa_pos", "alt_aa"])
    )

    # Plot correlation between time points
    dataset1_r1 = median_df.index[0]
    graphdf = median_df.loc[dataset1_r1].reset_index()[
        selcoeff_cols + ["confidence_score"]
    ]
    plot_timepoint_corr(graphdf, timepointsplot_outpath, plot_formats)

    median_long = median_df.melt(
        id_vars=mutation_attributes_aa,
        value_vars=selcoeff_cols,
        var_name="Compared timepoints",
        value_name="s",
        ignore_index=False,
    ).reset_index()

    # Output dataframe to plot more graphs (aggregating over all sample groups)
    median_long.to_csv(aa_df_outpath, index=False)

    # Calculate median across replicates for high confidence variants
    gby_score = median_df.groupby(["Replicate", "confidence_score"]).size()
    totseq = median_df.groupby("Replicate").size()
    perc_by_score = (gby_score / totseq).to_frame("proportion").reset_index()
    if (
        perc_by_score.loc[(perc_by_score.confidence_score == 1, "proportion")] < 0.75
    ).any():
        warnings.warn(
            f"Warning.. Group {sample_group} shows less than 75% of variants labeled with high confidence.\n"
            "Because only these variants are used to calculate a median score across replicates,"
            "you may want to double check the config file and adjust the rc_threshold parameter."
        )

    avg_df = (
        median_df[median_df.confidence_score == 1]
        .groupby(mutation_attributes_aa)[selcoeff_cols]
        .agg(
            [
                "median",
                lambda x: np.percentile(x, 2.5),
                lambda x: np.percentile(x, 97.5),
            ]
        )
    )

    # Rename columns
    cols_to_rename = [x for x in avg_df.columns if x not in mutation_attributes_aa]
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
    snakemake.config["samples"]["path"],
    snakemake.config["barcode"]["attributes"],
    snakemake.config["barcode"]["rc_level"],
    snakemake.config["filter"]["rc_threshold"],
    snakemake.config["rc_aims"]["exp_rc_per_var"],
    [x for x in snakemake.config["plots"]["format"] if x != "svg"],
)
