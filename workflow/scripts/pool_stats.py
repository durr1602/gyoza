#from snakemake.script import snakemake
from plotting_functions import plot_stacked_barplot, plot_unexp_plot
#from scripts.plotting_functions import plot_stacked_barplot, plot_unexp_plot
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
        list_df.append(pd.read_csv(f)[["Total_raw_reads", "Trimming", "Merging", "Aggregating"]])
    
    stats_df = pd.concat(list_df, ignore_index=True)
    
    list_df = []
    for f in sample_unexpected:
        # Retrieve sample name from file name
        sample_name = f.split('_unexpected.csv')[0]
        sample_unexp = pd.read_csv(f)[["readcount"]]
        sample_unexp["Sample_name"] = sample_name
        list_df.append(sample_unexp)

    unexp_all_seqs = pd.concat(list_df, ignore_index=True)

    plot_unexp_plot(unexp_all_seqs, unexpplot_outpath, plot_formats)

    unexp_df = (unexp_all_seqs
                .groupby("Sample_name")[["readcount"]]
                .sum()
                )

    stacked_df = pd.concat([stats_df,
                            unexp_df.rename(columns={"readcount": "Unexpected"})
                            ], axis=1
                            )

    stacked_df["OK"] = stacked_df["Total_raw_reads"] - stacked_df[
        ["Trimming", "Merging", "Aggregating", "Unexpected"]].sum(axis=1)

    stacked_df.sort_index(inplace=True)

    stacked_df.to_csv(csv_outpath)

    plot_stacked_barplot(stacked_df, barplot_outpath, exp_rc_per_sample, plot_formats)

    return

"""
get_pooled_stats(snakemake.input.read_stats,
                 snakemake.input.unexpected,
                 snakemake.output.all_stats,
                 snakemake.output.rc_filter_plot,
                 snakemake.output.unexp_rc_plot,
                 float(snakemake.config["rc_aims"]["exp_rc_per_sample"]),
                 [x for x in snakemake.config["plots"]["format"] if x != "svg"],
                 )
"""
import os
get_pooled_stats([f for f in os.listdir("/home/rodur28/git_repos/gyoza/results/stats/")],
                 [f for f in os.listdir("/home/rodur28/git_repos/gyoza/results/df/unexpected_seqs/")],
                  "/home/rodur28/git_repos/gyoza/results/test.csv",
                  "/home/rodur28/git_repos/gyoza/results/test1.svg",
                  "/home/rodur28/git_repos/gyoza/results/test2.svg",
                  3e4, ["svg"]
                 )