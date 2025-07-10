rule plot_heatmaps:
    input:
        "results/df/avg_scores_{group_key}.csv",
    output:
        report(
            "results/graphs/heatmap_fitness_{group_key}_{t}.svg",
            "../report/heatmap_aa.rst",
            category="3. Functional impact",
            subcategory="3.3. Heatmaps of functional impact",
            labels={"figure": "{group_key}_{t}"},
        ),
    params:
        wtaa=lambda w: extract_param_from_header(
            w, [f"results/df/avg_scores_{w.group_key}.csv"]
        ),
        level="aa",
    message:
        f"Plotting heatmaps of raw read counts.."
    log:
        "logs/8_scores/plot_heatmaps_{group_key}_{t}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/plot_heatmaps.py"
