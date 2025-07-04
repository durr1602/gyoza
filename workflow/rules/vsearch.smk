rule vsearch_fastx_uniques:
    input:
        fastx_uniques="results/2_merge/{sample}_merged.fasta",
    output:
        fastaout=temp("results/3_aggregate/{sample}_aggregated.fasta"),
    message:
        "Counting reads for every unique sequence..."
    log:
        "logs/3_aggregate/vsearch-sample={sample}.stats",
    params:
        # --minuniquesize discards sequences with abundance inferior to value
        # --relabel new headers with specified string and ticker (1,2,3...)
        # --sizeout conserve abundance annotations
        extra="--minuniquesize 2 --relabel seq --sizeout",
    envmodules:
        # If to be used, update the following, run module avail to see installed modules and versions
        "vsearch/2.29.3",
    wrapper:
        "v5.8.0/bio/vsearch"
