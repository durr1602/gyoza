rule annotate_random:
    input:
        rules.parse_fasta.output,
    output:
        annot_rc="results/df/unfiltered/{sample}_unfiltered-rc.csv",
        indels="results/df/indels/{sample}_indels.csv",
    params:
        position_offset=lambda wildcards: sample_layout.loc[
            wildcards.sample, "Pos_start"
        ],
    message:
        f"Annotating mutants observed in sequencing data"
    log:
        "logs/6_annotate/annotate_random-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/annotate_mutants.py"


rule filter_random_mutants:
    input:
        rules.annotate_random.output["annot_rc"],
    output:
        filtered="results/df/annotated_readcounts/{sample}_annot_rc.csv",
        discarded="results/df/unexpected_seqs/{sample}_unexpected.csv",
    message:
        f"Filtering out mutants based on number of amino acid changes"
    log:
        "logs/6_annotate/filter_random-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/filter_random.py"
