rule generate_report:
    input:
        "results/graphs/rc_filter_plot.svg",
        "results/graphs/unexp_rc_plot.svg",
        "results/graphs/rc_var_plot.svg",
        "results/graphs/scoeff_violin_plot.svg",
        "results/graphs/replicates_heatmap_plot.svg",
        "results/graphs/replicates_plot.svg",
        "results/graphs/s_through_time_plot.svg",
        expand("results/graphs/hist_plot_{group_key}.svg", group_key=REPORTED_GROUPS),
        expand("results/graphs/upset_plot_{group_key}.svg", group_key=REPORTED_GROUPS),
        expand(
            "results/graphs/timepoints_plot_{group_key}.svg", group_key=REPORTED_GROUPS
        ),
    output:
        "results/report.html",
    message:
        "Generate HTML report..."
    log:
        "logs/9_report/report.log"
    run:
        # Get list of successfully generated inputs
        completed = [f for f in input if Path(f).exists()]

        # Add QC report if applicable
        if config["qc"]["perform"] and Path("results/0_qc/multiqc.html").exists():
            completed.append("multiqc")

        if completed:
            shell(f"snakemake {' '.join(completed)} --report {output}")

        else:
            print(
                ">> No completed steps to include in report. Skipping report generation."
            )
