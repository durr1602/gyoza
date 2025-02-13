rule vsearch:
    input:
        rules.pandaseq.output
    output:
        'results/3_aggregate/{sample}_aggregated.fasta'
    resources:
        threads = 1, # This command of vsearch is not multi-threaded
        time = lambda _, input, attempt: max(0.0002*input.size_mb + (attempt-1)*0.0002*input.size_mb, 1)
    message:
        "Counting reads for every unique sequence..."
    log:
        'logs/3_aggregate/vsearch-sample={sample}.stats'
    conda:
        '../envs/vsearch.yaml'
    envmodules:
        # Update the following, run module avail to see installed modules and versions
        'vsearch/2.28.1'
    shell:
        ## Flags for vsearch
        # --fastx_uniques aggregates identical sequences
        # --minuniquesize discards sequences with abundance inferior to value
        # --relabel new headers with specified string and ticker (1,2,3...)
        # --sizeout conserve abundance annotations
        r"""
        vsearch --fastx_uniques {input} \
        --minuniquesize 2 \
        --relabel seq \
        --sizeout \
        --fastaout {output} &> {log}
        """