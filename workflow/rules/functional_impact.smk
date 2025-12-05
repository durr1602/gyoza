rule estimate_functional_impact:
    input:
        freq_df=rules.process_read_counts.output.freq_df,
        nbgen=NBGEN_PATH,
    output:
        selcoeffs=temp("results/df/all_scores_{group_key}.csv"),
        avg_scores=temp("results/df/avg_scores_{group_key}.csv"),
        timepoints_plot=report(
            "results/graphs/timepoints_plot_{group_key}.svg",
            "../report/timepoints_plot.rst",
            category="3. Functional impact",
            subcategory="3.1. Correlation between time points",
            labels={"figure": "{group_key}"},
        ),
        aa_df="results/df/agg_aa/aa_{group_key}.csv",
    params:
        all_attributes=SAMPLE_ATTR + SCREEN_ATTR,
        barcode_attributes=BC_ATTR,
        reported_groups=REPORTED_GROUPS_WITH_OUTPUTS,
        plot_formats=[x for x in config["plot_formats"] if x != "svg"],
    message:
        "Processing allele frequencies... converting to functional impact scores"
    log:
        "logs/9_scores/process_read_counts_{group_key}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/estimate_impact.py"


rule plot_scores:
    input:
        aa_df=expand(
            rules.estimate_functional_impact.output.aa_df,
            group_key=REPORTED_GROUPS_WITH_OUTPUTS,
        ),
    output:
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
        "logs/9_scores/plot_scores.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/plot_scores.py"


rule aggregate_dfs:
    input:
        all_df=expand(
            rules.estimate_functional_impact.output.selcoeffs,
            group_key=ATTR_GROUPS_WITH_OUTPUTS,
        ),
        avg_df=expand(
            rules.estimate_functional_impact.output.avg_scores,
            group_key=ATTR_GROUPS_WITH_OUTPUTS,
        ),
    output:
        selcoeffs="results/df/all_scores.csv",
        avg_scores="results/df/avg_scores.csv",
    message:
        "Exporting final dataframes"
    log:
        "logs/9_scores/aggregate_dfs.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/aggregate_dfs.py"
