rule stats:
    input:
        cutadapt_logs = expand(rules.cutadapt.log, sample=samples),
        pandaseq_logs = expand(rules.pandaseq.log, sample=samples),
        vsearch_logs = expand(rules.vsearch.log, sample=samples)
    output:
        'results/read_stats.csv'
    resources:
        threads = 1,
        time = lambda _, attempt: f'00:{attempt}:00'
    message:
        "Parsing log files to aggregate read statistics..."
    log:
        'logs/4_stats/stats.log'
    conda:
        '../envs/jupyter_basic.yaml'
    script:
        '../scripts/get_read_stats.py'
