rule pandaseq:
    input:
        read1 = rules.cutadapt.output.read1,
        read2 = rules.cutadapt.output.read2
    output:
        temp('results/2_merge/{sample}_merged.fasta')
    resources:
        threads = 4,
        time = lambda _, input, attempt: max(0.002*input.size_mb + (attempt-1)*0.002*input.size_mb, 1)
    message:
        "Merging reads for {input.read1} and {input.read2}"
    log:
        'logs/2_merge/pandaseq-sample={sample}.stats'
    conda:
        '../envs/pandaseq.yaml'
    envmodules:
        # Update the following, run module avail to see installed modules and versions
        'pandaseq/2.11'
    shell:
        ## Flags for pandaseq
        # -O max overlap, important,related to Aviti sequencing tech
        # -k number of k-mers
        # -B allow input sequences to lack a barcode/tag
        # -t minimum threshold for alignment score (0-1)
        # -T number of threads, important, see doc
        # -d flags to decide what to include in the output log, see doc
        # -w output
        r"""
        pandaseq -f {input.read1} -r {input.read2} \
        -O 625 \
        -k 4 \
        -B \
        -t 0.6 \
        -T {threads} \
        -d bFSrk \
        -w {output} &> {log}
        """
