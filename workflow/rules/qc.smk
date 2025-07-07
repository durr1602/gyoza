#### Adapted from https://github.com/akcorut/kGWASflow/

# =================================================================================================
#     Generate QC Stats Using FastQC on Raw Reads
# =================================================================================================


rule fastqc:
    input:
        lambda wildcards: f'{config["reads"]["path"]}{sample_layout.loc[wildcards.sample, wildcards.RF]}',
    output:
        html="results/0_qc/{sample}_{RF}_fastqc.html",
        zip="results/0_qc/{sample}_{RF}_fastqc.zip",
    log:
        "logs/0_qc/{sample}_{RF}_fastqc.log",
    message:
        "Performing quality control analysis using FastQC on the following file: {input}"
    wrapper:
        "v5.0.2/bio/fastqc"


# =================================================================================================
#     MultiQC
# =================================================================================================


rule multiqc:
    input:
        expand("results/0_qc/{sample}_{RF}_fastqc.zip", sample=SAMPLES, RF=RF_VALS),
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
    message:
        "Aggregating FastQC results with MultiQC..."
    wrapper:
        "v5.0.2/bio/multiqc"


# =================================================================================================
