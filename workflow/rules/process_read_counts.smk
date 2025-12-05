readcounts_by_group = {
    k: [f"results/df/annotated_readcounts/{s}_annot_rc.csv" for s in v]
    for k, v in final_groups_str.items()
}


rule process_read_counts:
    input:
        LAYOUT_PATH,
        readcounts=lambda wildcards: readcounts_by_group.get(wildcards.group_key, []),
    output:
        freq_df="results/df/distribution_freq/freq_{group_key}.csv",
        hist_plot=report(
            "results/graphs/hist_plot_{group_key}.svg",
            "../report/hist_plot.rst",
            category="2. Read counts and allele frequencies",
            subcategory="2.1. {group_key}",
            labels={"figure": "2.1.a. Raw read count per variant"},
        ),
        upset_plot=report(
            "results/graphs/upset_plot_{group_key}.svg",
            "../report/upset_plot.rst",
            category="2. Read counts and allele frequencies",
            subcategory="2.1. {group_key}",
            labels={"figure": "2.1.b. Overlap across time points and replicates"},
        ),
    params:
        layout=sample_layout,
        reported_samples=REPORTED_SAMPLES,
        readcount_level=RC_LEVEL,
        barcode_attributes=BC_ATTR,
        rc_threshold=config["reads"]["rc_threshold"],
        plot_formats=[x for x in config["plot_formats"] if x != "svg"],
    message:
        "Processing read counts... converting to allele frequencies"
    log:
        "logs/8_freqs/process_read_counts_{group_key}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/process_rc.py"


rule plot_freq:
    input:
        freq_df=expand(rules.process_read_counts.output.freq_df, group_key=ATTR_GROUPS),
    output:
        rc_var_plot=report(
            "results/graphs/rc_var_plot.svg",
            "../report/rc_var_plot.rst",
            category="2. Read counts and allele frequencies",
            subcategory="2.2. Aggregated",
            labels={"figure": "2.2.a. Distribution of allele frequencies"},
        ),
    params:
        reported_samples=REPORTED_SAMPLES,
        readcount_level=RC_LEVEL,
        plot_formats=[x for x in config["plot_formats"] if x != "svg"],
    message:
        "Plot distribution of allele frequencies"
    log:
        "logs/8_freqs/plot_freq.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/plot_freqs.py"
