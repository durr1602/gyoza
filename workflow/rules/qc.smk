rule fastp:
    input:
        sample=lambda w: [f"{READS_PATH}/{f}" for f in sample_layout.loc[w.sample, RF_VALS]]
    output:
        html="results/0_qc/{sample}.fastp.html",
        json="results/0_qc/{sample}.fastp.json",
    log:
        "logs/0_qc/{sample}_fastp.log",
    params:
        # --dont_eval_duplication rate (much faster, DMS data is duplicated data)
        extra="--dont_eval_duplication"
    message:
        "Performing quality control analysis using Fastp on the following file: {input}"
    wrapper:
        "v7.7.0/bio/fastp"


# Using an older MultiQC/wrapper version -> lighter env
rule multiqc:
    input:
        expand(
            "results/0_qc/{sample}.fastp.json",
            sample=REPORTED_SAMPLES,
        ),
    output:
        report(
            "results/0_qc/multiqc.html",
            "../report/qc.rst",
            category="0. Quality control",
            labels={"report": "Interactive QC report"},
        ),
    log:
        "logs/0_qc/multiqc.log",
    params:
        extra="-v -d --interactive",
        use_input_files_only=True,
    message:
        "Aggregating Fastp results with MultiQC..."
    wrapper:
        "v5.0.2/bio/multiqc"
