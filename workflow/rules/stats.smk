rule stats:
    input:
        cutadapt_logs = expand(rules.cutadapt.output.qc, sample=samples),
        pandaseq_logs = expand(rules.pandaseq.log, sample=samples),
        vsearch_logs = expand(rules.vsearch.log, sample=samples)
    output:
        'results/read_stats.csv'
    resources:
        threads = 1,
        time = lambda _, input, attempt: max((0.001*input.size_mb + (attempt-1)*0.001*input.size_mb).__ceil__(), 1)
    message:
        "Parsing log files to aggregate read statistics..."
    log:
        'logs/4_stats/stats.log'
    conda:
        '../envs/jupyter.yaml'
    script:
        '../scripts/get_read_stats.py'
