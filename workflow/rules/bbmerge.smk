rule bbmerge:
    input:
        in1=rules.cutadapt.output.fastq1,
        in2=rules.cutadapt.output.fastq2,
    output:
        out=temp("results/2_merge/{sample}_merged.fasta"),
    message:
        "Merging reads for {input.in1} and {input.in2}"
    log:
        "logs/2_merge/bbmerge-sample={sample}.stats",
    params:
        command="bbmerge.sh",
    threads: 8
    wrapper:
        "v8.0.3/bio/bbtools"  # bbtools v39.52
