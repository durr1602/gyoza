rule stats:
    input:
        cutadapt_log="logs/1_trim/cutadapt-sample={sample}.stats",
        pandaseq_log="logs/2_merge/pandaseq-sample={sample}.stats",
        vsearch_log="logs/3_aggregate/vsearch-sample={sample}.stats",
        N_discarded_log="logs/4_readcounts/parse-fasta-sample={sample}.log",
    output:
        temp("results/stats/stats-sample={sample}.csv"),
    params:
        is_paired=config["reads"]["paired"],
    message:
        "Parsing log files to aggregate read statistics..."
    log:
        "logs/7_stats/stats-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/get_read_stats.py"
