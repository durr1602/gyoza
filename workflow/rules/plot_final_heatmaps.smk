rule format_scores:
    input:
        "results/df/avg_scores_{group_key}.csv",
    output:
        heatmap_df="results/df/formatted_scores/{group_key}_{t}_format_s.csv",
        heatmap_meta="results/heatmap_metadata/{group_key}_{t}_s.pkl",
    params:
        position_offset=lambda wildcards: pos_offset_by_group[wildcards.group_key],
    message:
        f"Format functional impact scores.."
    log:
        "logs/9_heatmaps/format_s_{group_key}_{t}",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/format_s.py"


rule plot_heatmaps:
    input:
        heatmap_df=rules.format_scores.output.heatmap_df,
        heatmap_meta=rules.format_scores.output.heatmap_meta,
    output:
        report(
            "results/graphs/heatmap_fitness_{group_key}_{t}.svg",
            "../report/heatmap_aa.rst",
            category="3. Functional impact",
            subcategory="3.3. Heatmaps of functional impact",
            labels={"figure": "{group_key}_{t}"},
        ),
    params:
        plot_formats=[x for x in config["plot_formats"] if x != "svg"],
    message:
        f"Plotting heatmaps of raw read counts.."
    log:
        "logs/9_heatmaps/plot_heatmaps_{group_key}_{t}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/plot_heatmaps.py"
