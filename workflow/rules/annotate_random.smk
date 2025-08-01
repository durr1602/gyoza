rule inject_WT:
    input:
        rules.parse_fasta.output,
    output:
        temp("results/df/readcounts/{sample}_rc_with_WT.csv"),
    params:
        wt=lambda wildcards: mutseq_to_wtseq[sample_to_mutseq[wildcards.sample]],
    log:
        "logs/6_annotate/inject_WT-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/inject_WT.py"


rule annotate_random:
    input:
        rules.inject_WT.output,
    output:
        annot_rc="results/df/unfiltered/{sample}_unfiltered-rc.csv",
        indels="results/df/indels/{sample}_indels.csv",
    params:
        position_offset=lambda wildcards: sample_layout.loc[
            wildcards.sample, "Pos_start"
        ],
        genetic_code=GEN_CODE_PATH,
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
    params:
        Nham_aa_max=int(config["random"]["Nham_aa_max"]),
    message:
        f"Filtering out mutants based on number of amino acid changes"
    log:
        "logs/6_annotate/filter_random-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/filter_random.py"
