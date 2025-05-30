rule stats:
    input:
        cutadapt_log=rules.cutadapt.output.qc,
        pandaseq_log=rules.pandaseq.log,
        vsearch_log=rules.vsearch_fastx_uniques.log,
    output:
        "results/stats/stats-sample={sample}.csv",
    message:
        "Parsing log files to aggregate read statistics..."
    log:
        "logs/7_stats/stats-sample={sample}.log",
    conda:
        "../envs/jupyter.yaml"
    script:
        "../scripts/get_read_stats.py"
