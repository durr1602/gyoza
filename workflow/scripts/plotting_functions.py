import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams["svg.fonttype"] = "none"

def plot_stacked_barplot(df, outpath, exp_rc_per_sample, plot_formats):
    samples = df.index.to_list()
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