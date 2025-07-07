rule pandaseq:
    input:
        read1=rules.cutadapt.output.fastq1,
        read2=rules.cutadapt.output.fastq2,
    output:
        temp("results/2_merge/{sample}_merged.fasta"),
    message:
        "Merging reads for {input.read1} and {input.read2}"
    log:
        "logs/2_merge/pandaseq-sample={sample}.stats",
    conda:
        "../envs/pandaseq.yaml"
    envmodules:
        # If to be used, update the following, run module avail to see installed modules and versions
        "pandaseq/2.11",
    shell:
        ## Flags for pandaseq
        # -O max overlap, important, related to Aviti sequencing tech
        # -k number of k-mers
        # -N exclude sequences N(s)
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
        -N \
        -t 0.6 \
        -T {threads} \
        -d bfsrk \
        -w {output} &> {log}
        """
