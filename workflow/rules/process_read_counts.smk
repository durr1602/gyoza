readcounts_by_group = {
    k: [f"results/df/annotated_readcounts/{s}_annot_rc.csv" for s in v]
    for k, v in grouped_samples_str.items()
}


rule process_read_counts:
    input:
        readcounts=lambda wildcards: readcounts_by_group.get(wildcards.group_key, []),
        nbgen=config["samples"]["generations"],
    output:
        selcoeffs="results/df/all_scores_{group_key}.csv",
        avg_scores="results/df/avg_scores_{group_key}.csv",
        wt_scores="results/df/avg_WT_{group_key}.csv",
        hist_plot=report(
            "results/graphs/hist_plot_{group_key}.svg",
            "../report/hist_plot.rst",
            category="2. Read processing",
            subcategory="2.1. {group_key}",
            labels={"figure": "2.1.a. Raw read count per variant"},
        ),
        upset_plot=report(
            "results/graphs/upset_plot_{group_key}.svg",
            "../report/upset_plot.rst",
            category="2. Read processing",
            subcategory="2.1. {group_key}",
            labels={"figure": "2.1.b. Overlap across time points and replicates"},
        ),
        timepoints_plot=report(
            "results/graphs/timepoints_plot_{group_key}.svg",
            "../report/timepoints_plot.rst",
            category="3. Functional impact",
            subcategory="3.2. Correlation between time points",
            labels={"figure": "{group_key}"},
        ),
        freq_df="results/df/distribution_freq/freq_{group_key}.csv",
        aa_df="results/df/agg_aa/aa_{group_key}.csv",
    message:
        "Processing read counts... converting to functional impact scores"
    log:
        "logs/8_scores/process_read_counts_{group_key}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/process_rc.py"


rule plot_scores:
    input:
        freq_df=expand(rules.process_read_counts.output.freq_df, group_key=GROUP_KEYS),
        aa_df=expand(rules.process_read_counts.output.aa_df, group_key=GROUP_KEYS),
    output:
        rc_var_plot=report(
            "results/graphs/rc_var_plot.svg",
            "../report/rc_var_plot.rst",
            category="2. Read processing",
            subcategory="2.2. Aggregated",
            labels={"figure": "2.2.a. Distribution of allele frequencies"},
        ),
        scoeff_violin_plot=report(
            "results/graphs/scoeff_violin_plot.svg",
            "../report/scoeff_violin_plot.rst",
            category="3. Functional impact",
            subcategory="3.1. Aggregated",
            labels={"figure": "3.1.a. Distribution of functional impact scores"},
        ),
        replicates_heatmap_plot=report(
            "results/graphs/replicates_heatmap_plot.svg",
            "../report/replicates_heatmap_plot.rst",
            category="3. Functional impact",
            subcategory="3.1. Aggregated",
            labels={"figure": "3.1.b. Correlation between replicates (1/2)"},
        ),
        replicates_plot=report(
            "results/graphs/replicates_plot.svg",
            "../report/replicates_plot.rst",
            category="3. Functional impact",
            subcategory="3.1. Aggregated",
            labels={"figure": "3.1.c. Correlation between replicates (2/2)"},
        ),
        s_through_time_plot=report(
            "results/graphs/s_through_time_plot.svg",
            "../report/s_through_time_plot.rst",
            category="3. Functional impact",
            subcategory="3.1. Aggregated",
            labels={"figure": "3.1.d. Functional impact over time"},
        ),
    message:
        "Aggregating dataframes to plot functional impact scores"
    log:
        "logs/8_scores/plot_scores.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/plot_scores.py"
