rule format_read_counts:
    input:
        "results/df/annotated_readcounts/{sample}_annot_rc.csv",
    output:
        heatmap_df="results/df/formatted_readcounts/{sample}_format_rc.csv",
        heatmap_meta="results/heatmap_metadata/{sample}_rc.pkl",
    params:
        exp_rc=float(config["rc_aims"]["exp_rc_per_sample"]),
    message:
        f"Format read counts.."
    log:
        "logs/6_format/format_rc-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/format_rc.py"


rule plot_read_counts:
    input:
        heatmap_df=rules.format_read_counts.output.heatmap_df,
        heatmap_meta=rules.format_read_counts.output.heatmap_meta,
    output:
        report(
            "results/graphs/heatmap_readcount_{sample}.svg",
            "../report/heatmap_rc.rst",
            category="1. Read filtering",
            subcategory="1.2. Heatmaps of raw read counts",
            labels={"figure": "{sample}"},
        ),
    message:
        f"Plotting heatmaps of raw read counts.."
    log:
        "logs/7_plot/plot_rc-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/plot_heatmaps.py"
