from snakemake.script import snakemake
from scripts.plotting_functions import plot_stacked_barplot, plot_unexp_plot
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"

def get_pooled_stats(sample_stats, sample_unexpected, csv_outpath, barplot_outpath, unexpplot_outpath, exp_rc_per_sample, plot_formats):
    """
    Pool read count statistics from all samples.
    Merges info with unexpected sequences (seen in sequencing dataset but not expected or filtered out)
    """
    list_df = []
    for f in sample_stats:
        list_df.append(pd.read_csv(f)[["Sample_name", "Total_raw_reads", "Trimming", "Merging", "Aggregating"]])
    
    stats_df = pd.concat(list_df, ignore_index=True)
    
    list_df = []
    for f in sample_unexpected:
        sample_unexp = pd.read_csv(f)[["Sample_name", "readcount"]]
        list_df.append(sample_unexp)

    unexp_all_seqs = pd.concat(list_df, ignore_index=True)

    plot_unexp_plot(unexp_all_seqs, unexpplot_outpath, plot_formats)

    unexp_df = (unexp_all_seqs
                .groupby("Sample_name")[["readcount"]]
                .sum()
                )

    stacked_df = pd.merge(left=stats_df,
                          right=unexp_df.rename(columns={"readcount": "Unexpected"}),
                          on='Sample_name'
                          )

    stacked_df["OK"] = stacked_df["Total_raw_reads"] - stacked_df[
        ["Trimming", "Merging", "Aggregating", "Unexpected"]].sum(axis=1)

    stacked_df.sort_index(inplace=True)

    stacked_df.to_csv(csv_outpath)

    plot_stacked_barplot(stacked_df, barplot_outpath, exp_rc_per_sample, plot_formats)

    return

get_pooled_stats(snakemake.input.read_stats,
                 snakemake.input.unexpected,
                 snakemake.output.all_stats,
                 snakemake.output.rc_filter_plot,
                 snakemake.output.unexp_rc_plot,
                 float(snakemake.config["rc_aims"]["exp_rc_per_sample"]),
                 [x for x in snakemake.config["plots"]["format"] if x != "svg"],
                 )