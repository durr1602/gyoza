rule annotate_mutants:
    input:
        rules.compare_to_sequencing.output.observed,
        LAYOUT_PATH,
        GEN_CODE_PATH,
    output:
        annot_rc="results/df/annotated_readcounts/{sample}_annot_rc.csv",
    params:
        position_offset=lambda wildcards: sample_layout.loc[
            wildcards.sample, "Pos_start"
        ],
        genetic_code=GEN_CODE,
    message:
        f"Annotating expected mutants (number of mutations, alternative codons, etc)"
    log:
        "logs/6_annotate/annotate_mutants-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/annotate_mutants.py"
