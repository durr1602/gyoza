rule cutadapt:
    input:
        lambda wildcards: [
            f"{config['reads']['path']}{sample_layout.loc[wildcards.sample, 'R1']}",
            f"{config['reads']['path']}{sample_layout.loc[wildcards.sample, 'R2']}",
        ],
    params:
        adapters=lambda wildcards: f"-g {sample_layout.loc[wildcards.sample, 'N_forward']} -G {sample_layout.loc[wildcards.sample, 'N_reverse']}",
        extra="-e 0.15 --no-indels --discard-untrimmed",
    output:
        fastq1=temp("results/1_trim/{sample}_trimmed.R1.fastq.gz"),
        fastq2=temp("results/1_trim/{sample}_trimmed.R2.fastq.gz"),
        qc="logs/1_trim/cutadapt-sample={sample}.stats",
    resources:
        threads=10,
        time=lambda _, input, attempt: max(
            (0.005 * input.size_mb + (attempt - 1) * 0.005 * input.size_mb).__ceil__(),
            1,
        ),
    message:
        "Trimming constant sequences from input files {input[0]} and {input[1]}"
    log:
        "logs/1_trim/cutadapt-sample={sample}.err",
    envmodules:
        # If to be used, update the following, run module avail to see installed modules and versions
        "cutadapt/5.0",
    wrapper:
        "v5.8.0/bio/cutadapt/pe"
