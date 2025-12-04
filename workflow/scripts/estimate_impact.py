"""Module to convert allele frequencies into functional impact scores."""

from snakemake.script import snakemake
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib

matplotlib.use("Agg")  # Non-GUI backend for generating plots without display
import matplotlib.pyplot as plt

plt.rcParams["svg.fonttype"] = "none"

MUTATION_ATTRIBUTES = [
    "mutated_codon",
    "mutation_aa_pos",
    "mutation_alt_codons",
    "mutation_alt_aa",
    "mutation_type",
]

PROT_SEQ_ATTRIBUTES = [
    "Nham_aa",
    "aa_seq",
]

MISSENSE_AA_ATTRIBUTES = [
    "aa_pos",
    "alt_aa",
    "wt_aa",
]

SEQUENCE_ATTRIBUTES = (
    [
        "nt_seq",
        "Nham_nt",
        "Nham_codons",
        "confidence_score",
    ]
    + PROT_SEQ_ATTRIBUTES
    + MISSENSE_AA_ATTRIBUTES
)


def aggregate_multiple_attr(g):
    r"""Returns not-applicable if more than 1 value, else returns value.

    Parameters
    ----------
    g : pandas.Series
        e.g. missense_aa_attributes

    Returns
    -------
    str
        Either "not-applicable" or unique value from ``g``
    """
    unique = g.unique()
    if len(unique) > 1:
        return "not-applicable"
    else:
        return unique[0]


def plot_timepoint_corr(df, outpath, plot_formats):
    r"""Plot pairwise comparisons of functional impact scores between time points.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe of functional impact scores.
        Should contain column ``Replicate``.
    outpath : str
        Path to save plot as SVG (should end with ``.svg``).
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.
    """
    # Check number of columns
    if len([x for x in df.columns if x not in ["Sample attributes", "Replicate"]]) <= 1:
        f, ax = plt.subplots(figsize=(4, 4))
        ax.text(0.5, 0.5, "Not enough time points to plot", ha="center", va="center")
        ax.set_axis_off()  # hide axes
    else:
        sample_group = df["Sample attributes"].values[0]
        g = sns.pairplot(
            df,
            hue="Replicate",
            palette="hls",
            plot_kws={"s": 8, "alpha": 0.2},
            height=1.5,
            corner=True,
        )
        g.tight_layout()
        plt.subplots_adjust(top=0.9)
    plt.suptitle(f"{sample_group}")

    plt.savefig(outpath, format="svg", dpi=300)
    [
        plt.savefig(f"{outpath.split('.svg')[0]}.{x}", format=x, dpi=300)
        for x in plot_formats
    ]
    return


def get_selcoeffs(
    data_path,
    nbgen_path,
    outpath,
    avg_outpath,
    timepointsplot_outpath,
    aa_df_outpath,
    all_attributes,
    barcode_attributes,
    plot_formats,
):
    r"""Convert allele frequencies into functional impact scores for matched samples.

    Parameters
    ----------
    data_path : str
        Path to CSV-formatted dataframe of allele frequencies.
        Should contain columns:

        * ``Sample attributes`` (**str**, serialized unique combinations of
        sample and screening attributes)
        * ``Timepoint`` (**str**)
        * ``Replicate`` (**str**)
        * ``frequency`` (**float**, allele frequency of the variant)

        All mutation attributes:

        * ``mutated_codon`` (**str**, ``1`` for the first mutated codon,
        ``2`` for the second mutated codon in the same ``nt_seq``, etc.)
        * ``mutation_aa_pos`` (**str**, position in the protein sequence at
        which the wild-type codon has been replaced by ``mutation_alt_codons``
        * ``mutation_alt_codons`` (**str**, mutated codon at position
        ``mutation_aa_pos``)
        * ``mutation_alt_aa`` (**str**, residue translated from
        ``mutation_alt_codons``)
        * ``mutation_type`` (**str**, either ``wt`` (no mutation), ``silent``
        (``mutation_alt_aa`` corresponds to wild-type residue), ``nonsense``
        (``mutation_alt_codons`` is a stop codon) or ``missense``)

        All sequence attributes:

        * ``nt_seq`` (**str**, nucleotide sequence of the variant)
        * ``confidence_score`` ({1, 2, 3}, ``1==high``, ``3==low``)
        * ``aa_seq`` (**str**, amino acid sequence of the variant)
        * ``Nham_nt`` (**str**, number of nucleotide changes when comparing
        ``nt_seq`` to the wild-type sequence)
        * ``Nham_aa`` (**str**, number of amino acid changes when comparing
        ``aa_seq`` to the wild-type sequence)
        * ``Nham_codons`` (**str**, number of codon changes when comparing
        ``nt_seq`` to the wild-type sequence)
        * ``aa_pos`` (**list** of positions in the protein sequence
        at which the wild-type residue has been mutated
        * ``alt_aa`` (**list** of alternative amino acid residues at ``aa_pos``)
        * ``wt_aa`` (**list** of wild-type amino acid residues at ``aa_pos``)

    nbgen_path : str
        Path to CSV-formatted dataframe containing the number of mitotic
        generations between T0 and each time point.
    outpath : str
        Path to save output dataframe of functional impact scores.
    avg_outpath : str
        Path to save output dataframe of fitness and error values (functional
        impact scores averaged over replicates for high confidence variants only).
    timepointsplot_outpath : str
        Path to save plot with comparisons between time points as SVG
        (should end with ``.svg``).
    aa_df_outpath : str
        Path to save output dataframe of functional impact scores aggregated at
        the protein level.
    all_attributes : list of str
        List of sample and screening attributes.
        The corresponding values should feature in column ``Sample_attributes``
        of `data`.
    barcode_attributes : list of str
        List of barcode attributes (includes `rc_level`).
    plot_formats : list of str
        Formats other than SVG in which the plot should be saved.

    Raises
    ------
    Exception
        In case of null sample depth.

    Notes
    -----
    Functional impact scores are obtained with a log ratio method:

    .. math:: s_v=\ \log_2{\left(\frac{c_{v,output}}{\sum\nolimits_{i} c_{i,\ output}}\right)}\ -\log_2{\left(\frac{c_{v,input}}{\sum\nolimits_{i} c_{i,\ input}}\right)}

    with :math:`c_v` being the raw read count of a variant + 1,
    "input" being T0 and "output" designating any post-screening time point.
    """
    # Import data
    longfreq = pd.read_csv(data_path)

    # Retrieve sample group (both serialized and deserialized)
    sample_group = longfreq["Sample attributes"].values[0]
    sample_group_tuple = tuple(sample_group.split("__"))

    # Retrieve list of time points
    timepoints = sorted(longfreq.Timepoint.unique())

    # Calculate Log2(fold-change) for every time point relative to T0
    freq_wide = longfreq.pivot(
        index=["Sample attributes"]
        + MUTATION_ATTRIBUTES
        + SEQUENCE_ATTRIBUTES
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
        nbgen_df[all_attributes].apply(tuple, axis=1) == sample_group_tuple
    ]
    nbgen_wide = nbgen_group.pivot(
        index="Replicate", columns="Timepoint", values="Nb_gen"
    )
    nbgen_wide.columns = [f"{x}_gen" for x in nbgen_wide.columns]
    gen_cols = nbgen_wide.columns
    lfc_wide = freq_wide.reset_index().merge(
        right=nbgen_wide.reset_index(), on="Replicate"
    )

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

    # Save metadata in df for simplicity
    s_wide[all_attributes] = sample_group_tuple

    # Export full dataframe
    s_wide[
        all_attributes
        + ["Replicate"]
        + SEQUENCE_ATTRIBUTES
        + MUTATION_ATTRIBUTES
        + barcode_attributes
        + selcoeff_cols
    ].to_csv(outpath, index=False)

    # Calculate median functional impact score (over synonymous codons),
    # for each replicate separately,
    # from high confidence variants ONLY
    scoeff_agg = dict(zip(selcoeff_cols, ["median"] * len(selcoeff_cols)))
    missense_agg = dict(
        zip(
            MISSENSE_AA_ATTRIBUTES,
            [aggregate_multiple_attr] * len(MISSENSE_AA_ATTRIBUTES),
        )
    )
    agg_dict = {**scoeff_agg, **missense_agg}
    median_df = (
        s_wide[s_wide.confidence_score == 1]
        .groupby(["Sample attributes", "Replicate"] + PROT_SEQ_ATTRIBUTES)[
            selcoeff_cols + MISSENSE_AA_ATTRIBUTES
        ]
        .agg(agg_dict)
        .reset_index(level=PROT_SEQ_ATTRIBUTES + ["Sample attributes"])
    )

    # Plot correlation between time points
    plot_timepoint_corr(
        median_df.reset_index()[["Sample attributes", "Replicate"] + selcoeff_cols],
        timepointsplot_outpath,
        plot_formats,
    )

    # Reshape
    median_long = median_df.melt(
        id_vars=["Sample attributes"] + PROT_SEQ_ATTRIBUTES + MISSENSE_AA_ATTRIBUTES,
        value_vars=selcoeff_cols,
        var_name="Compared timepoints",
        value_name="s",
        ignore_index=False,
    ).reset_index()
    # Rename column to keep only output time point (all are compared relative to T0)
    median_long["Compared timepoints"] = median_long["Compared timepoints"].apply(
        lambda x: x.split("_")[1]
    )

    # Output dataframe to plot more graphs (aggregating over sample groups)
    median_long.to_csv(aa_df_outpath, index=False)

    # Calculate median across replicates for high confidence variants
    avg_df = median_df.groupby(PROT_SEQ_ATTRIBUTES + MISSENSE_AA_ATTRIBUTES)[
        selcoeff_cols
    ].agg(
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

    # Rename columns
    cols_to_rename = [
        x
        for x in avg_df.columns
        if x not in PROT_SEQ_ATTRIBUTES + MISSENSE_AA_ATTRIBUTES
    ]
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

    # Save metadata in df
    avg_df[all_attributes] = sample_group_tuple

    # Export dataframe with fitness and error values
    avg_df.reset_index()[
        all_attributes + PROT_SEQ_ATTRIBUTES + MISSENSE_AA_ATTRIBUTES + new_names
    ].to_csv(avg_outpath, index=False)

    return


get_selcoeffs(
    snakemake.input.freq_df,
    snakemake.input.nbgen,
    snakemake.output.selcoeffs,
    snakemake.output.avg_scores,
    snakemake.output.timepoints_plot,
    snakemake.output.aa_df,
    snakemake.params.all_attributes,
    snakemake.params.barcode_attributes,
    snakemake.params.plot_formats,
)
