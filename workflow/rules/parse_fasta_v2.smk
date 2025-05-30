rule parse_fasta:
    input:
        rules.vsearch_fastx_uniques.output,
    output:
        "results/df/readcounts/{sample}_rc.csv",
    params:
        lambda wildcards: mutseq_to_wtseq[sample_to_mutseq[wildcards.sample]]
    message:
        "Parsing fasta files to get read counts..."
    log:
        "logs/4_readcounts/parse-fasta-sample={sample}.log",
    conda:
        "../envs/jupyter.yaml"
    script:
        "../scripts/parse_fasta.py"
