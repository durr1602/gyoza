from snakemake.script import snakemake
from scripts.my_functions import concatenate_df
from scripts.plotting_functions import (
    plot_allele_freq,
    plot_scoeff_violin,
    plot_impact_over_time,
    plot_spearman_heatmaps,
    plot_replicate_scatter,
)
import pandas as pd

mutation_attributes_aa = [
    "mutated_codon",
    "aa_pos",
    "alt_aa",
    "Nham_aa",
    "mutation_type",
]


def get_allele_freq_plot(df_files, outpath, plot_formats):
    df = concatenate_df(df_files)
    mean_exp_freq = (
        df.groupby("Sample attributes")[["Mean_exp_freq"]].first().mean(axis=None)
    )
    plot_allele_freq(df, outpath, mean_exp_freq, plot_formats)
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
        index=mutation_attributes_aa + ["Sample attributes", "Compared timepoints"],
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
