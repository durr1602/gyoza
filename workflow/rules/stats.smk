rule stats:
    input:
        cutadapt_logs=expand(rules.cutadapt.output.qc, sample=samples),
        pandaseq_logs=expand(rules.pandaseq.log, sample=samples),
        vsearch_logs=expand(rules.vsearch_fastx_uniques.log, sample=samples),
    output:
        "results/read_stats.csv",
    message:
        "Parsing log files to aggregate read statistics..."
    log:
        "logs/4_stats/stats.log",
    conda:
        "../envs/jupyter.yaml"
    script:
        "../scripts/get_read_stats.py"
