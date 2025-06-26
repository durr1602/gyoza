report: "report/workflow.rst"


EXPECTED_GRAPHS = (
    [
        Path("results/graphs/rc_filter_plot.svg"),
        Path("results/graphs/unexp_rc_plot.svg"),
        Path("results/graphs/rc_var_plot.svg"),
        Path("results/graphs/scoeff_violin_plot.svg"),
        Path("results/graphs/replicates_heatmap_plot.svg"),
        Path("results/graphs/replicates_plot.svg"),
        Path("results/graphs/s_through_time_plot.svg"),
    ]
    + [Path(f"results/graphs/hist_plot_{k}.svg") for k in REPORTED_GROUPS]
    + [Path(f"results/graphs/upset_plot_{k}.svg") for k in REPORTED_GROUPS]
    + [Path(f"results/graphs/timepoints_plot_{k}.svg") for k in REPORTED_GROUPS]
)


checkpoint discover_graphs:
    """
    Checkpoint to detect which graphs have been generated.
    Ensures report gets generated even upon intermediate workflow run.
    """
    output:
        "results/graphs/.check"
    run:
        Path("results/graphs").mkdir(parents=True, exist_ok=True)
        Path(output[0]).touch()


def collect_graphs(wildcards):
    graph_dir = Path("results/graphs")
    collected = [f for f in EXPECTED_GRAPHS if (graph_dir / f.name).exists()]
    return collected


rule generate_report:
    input:
        collect_graphs,
    output:
        "results/report.html",
    message:
        "Generate HTML report..."
    log:
        "logs/9_report/report.log",
    run:
        # Get list of successfully generated inputs
        completed = list(input)

        # Add QC report if applicable
        qc_path = Path("results/0_qc/multiqc.html")
        if config["qc"]["perform"] and qc_path.exists():
            completed.append(qc_path)

        if completed:
            shell(f"snakemake {' '.join(str(f) for f in completed)} --report {output}")
            Path(str(output[0])).touch()  # unfortunately essential if we want a DAG-included dynamic report

        else:
            print(">> No expected graph found. Skipping report generation.")
