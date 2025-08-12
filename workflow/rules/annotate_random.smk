rule discard_indels:
    input:
        rules.parse_fasta.output.readcounts,
    output:
        observed="results/df/observed_seqs/{sample}_observed.csv",
        indels=temp("results/df/indels/{sample}_indels.csv"),
    params:
        wt=lambda wildcards: mutseq_to_wtseq[sample_to_mutseq[wildcards.sample]],
    log:
        "logs/6_annotate/discard_indels-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/discard_indels.py"


rule annotate_random:
    input:
        rules.discard_indels.output.observed,
    output:
        annot_rc="results/df/unfiltered/{sample}_unfiltered-rc.csv",
    params:
        position_offset=lambda wildcards: sample_layout.loc[
            wildcards.sample, "Pos_start"
        ],
        genetic_code=GEN_CODE_PATH,
    message:
        f"Annotating mutants observed in sequencing data (number of mutations, alternative codons, etc)"
    log:
        "logs/6_annotate/annotate_random-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/annotate_mutants.py"


rule filter_random_mutants:
    input:
        no_indels=rules.annotate_random.output.annot_rc,
        indels=rules.discard_indels.output.indels,
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
