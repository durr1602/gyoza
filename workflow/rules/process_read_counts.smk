readcounts_by_group = {
    k: [f"results/df/annotated_readcounts/{s}_annot_rc.csv" for s in v]
    for k, v in final_groups_str.items()
}


rule process_read_counts:
    input:
        readcounts=lambda wildcards: readcounts_by_group.get(wildcards.group_key, []),
        nbgen=NBGEN_PATH,
    output:
        selcoeffs=temp("results/df/all_scores_{group_key}.csv"),
        avg_scores=temp("results/df/avg_scores_{group_key}.csv"),
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
            subcategory="3.1. Correlation between time points",
            labels={"figure": "{group_key}"},
        ),
        freq_df="results/df/distribution_freq/freq_{group_key}.csv",
        aa_df="results/df/agg_aa/aa_{group_key}.csv",
    params:
        layout=LAYOUT_PATH,
        sample_attributes=SAMPLE_ATTR,
        readcount_level=RC_LEVEL,
        barcode_attributes=BC_ATTR,
        rc_threshold=config["reads"]["rc_threshold"],
        plot_formats=[x for x in config["plot_formats"] if x != "svg"],
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
        freq_df=expand(
            rules.process_read_counts.output.freq_df,
            group_key=REPORTED_GROUPS_WITH_OUTPUTS,
        ),
        aa_df=expand(
            rules.process_read_counts.output.aa_df,
            group_key=REPORTED_GROUPS_WITH_OUTPUTS,
        ),
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
            subcategory="3.2. Aggregated",
            labels={"figure": "3.2.a. Distribution of functional impact scores"},
        ),
        replicates_heatmap_plot=report(
            "results/graphs/replicates_heatmap_plot.svg",
            "../report/replicates_heatmap_plot.rst",
            category="3. Functional impact",
            subcategory="3.2. Aggregated",
            labels={"figure": "3.2.b. Correlation between replicates (1/2)"},
        ),
        replicates_plot=report(
            "results/graphs/replicates_plot.svg",
            "../report/replicates_plot.rst",
            category="3. Functional impact",
            subcategory="3.2. Aggregated",
            labels={"figure": "3.2.c. Correlation between replicates (2/2)"},
        ),
        s_through_time_plot=report(
            "results/graphs/s_through_time_plot.svg",
            "../report/s_through_time_plot.rst",
            category="3. Functional impact",
            subcategory="3.2. Aggregated",
            labels={"figure": "3.2.d. Functional impact over time"},
        ),
    params:
        plot_formats=[x for x in config["plot_formats"] if x != "svg"],
    message:
        "Aggregating dataframes to plot functional impact scores"
    log:
        "logs/8_scores/plot_scores.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/plot_scores.py"


rule aggregate_dfs:
    input:
        all_df=expand(
            rules.process_read_counts.output.selcoeffs,
            group_key=ATTR_GROUPS_WITH_OUTPUTS,
        ),
        avg_df=expand(
            rules.process_read_counts.output.avg_scores,
            group_key=ATTR_GROUPS_WITH_OUTPUTS,
        ),
    output:
        selcoeffs="results/df/all_scores.csv",
        avg_scores="results/df/avg_scores.csv",
    message:
        "Exporting final dataframes"
    log:
        "logs/8_scores/aggregate_dfs.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/aggregate_dfs.py"
