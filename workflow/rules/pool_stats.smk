rule pool_stats:
    input:
        read_stats=expand(rules.stats.output, sample=SAMPLES),
        unexpected=expand(
            "results/df/unexpected_seqs/{sample}_unexpected.csv", sample=SAMPLES
        ),
    output:
        all_stats="results/all_stats.csv",
        rc_filter_plot=report(
            "results/graphs/rc_filter_plot.svg",
            "../report/rc_filter_plot.rst",
            category="1. Read filtering",
            labels={"figure": "1.1. Summary of filtered reads"},
        ),
        unexp_rc_plot=report(
            "results/graphs/unexp_rc_plot.svg",
            "../report/unexp_rc_plot.rst",
            category="1. Read filtering",
            labels={"figure": "1.2. Read counts of unexpected variants"},
        ),
    message:
        "Pooling statistics from all samples..."
    log:
        "logs/7_stats/pool_stats.log",
    conda:
        "../envs/jupyter.yaml"
    script:
        "../scripts/pool_stats.py"
