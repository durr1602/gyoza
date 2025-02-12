rule cutadapt:
    input:
        read1 = lambda wildcards: f"{config["reads"]["path"]}{sample_layout.loc[wildcards.sample, 'R1']}",
        read2 = lambda wildcards: f"{config["reads"]["path"]}{sample_layout.loc[wildcards.sample, 'R2']}"
    params:
        Nfwd = lambda wildcards: sample_layout.loc[wildcards.sample, "N_forward"],
        Nrev = lambda wildcards: sample_layout.loc[wildcards.sample, "N_reverse"]
    output:
        read1 = temp('results/1_trim/{sample}_trimmed.R1.fastq.gz'),
        read2 = temp('results/1_trim/{sample}_trimmed.R2.fastq.gz')
    resources:
        threads = 10,
        time = lambda _, input, attempt: max(0,5*input.size_mb + (attempt-1), attempt)
    message:
        "Trimming constant sequences forward={params.Nfwd} and reverse={params.Nrev} from input files {input.read1} and {input.read2}, respectively"
    log:
        'logs/1_trim/cutadapt-sample={sample}.stats'
    conda:
        '../envs/cutadapt.yaml'
    envmodules:
        # Update the following, run module avail to see installed modules and versions
        'cutadapt/4.9'
    shell:
        r"""
        cutadapt -e 0.15 \
        --no-indels \
        --discard-untrimmed \
        --cores={threads} \
        -g {params.Nfwd} -G {params.Nrev} \
        -o {output.read1} -p {output.read2} \
        {input.read1} {input.read2} &> {log}
        """
