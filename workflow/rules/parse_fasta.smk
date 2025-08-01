rule parse_fasta:
    input:
        rules.vsearch_fastx_uniques.output,
    output:
        "results/df/readcounts/{sample}_rc.csv",
    params:
        mutseq=lambda wildcards: sample_to_mutseq[wildcards.sample],
        readcount_level=RC_LEVEL,
    message:
        "Parsing fasta files to get read counts..."
    log:
        "logs/4_readcounts/parse-fasta-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/parse_fasta.py"
