rule cutadapt_se:
    input:
        lambda wildcards: f"{READS_PATH}/{sample_layout.loc[wildcards.sample, 'R1']}",
    params:
        adapters=lambda wildcards: f"-g {sample_layout.loc[wildcards.sample, 'N_forward']}",
        extra="-e 0.15 --no-indels --discard-untrimmed",
    output:
        touch("logs/2_merge/pandaseq-sample={sample}.stats"),
        fastq=temp("results/2_merge/{sample}_merged.fasta"),
        qc="logs/1_trim/cutadapt-sample={sample}.stats",
    message:
        "Trimming constant sequences from input file {input[0]}"
    log:
        "logs/1_trim/cutadapt-sample={sample}.err",
    envmodules:
        # If to be used, update the following, run module avail to see installed modules and versions
        "cutadapt/5.0",
    wrapper:
        "v7.1.0/bio/cutadapt/se"
