rule annotate_random:
    input:
        rules.parse_fasta.output,
    output:
        annot_rc="results/df/unfiltered/{sample}_unfiltered-rc.csv",
        indels="results/df/indels/{sample}_indels.csv",
    message:
        f"Annotating mutants observed in sequencing data"
    log:
        "logs/6_annotate/annotate_random-sample={sample}.log",
    conda:
        "../envs/jupyter.yaml"
    script:
        "../scripts/annotate_mutants.py"


rule filter_random_mutants:
    input:
        rules.annotate_random.output["annot_rc"],
    output:
        filtered="results/df/annotated_readcounts/{sample}_annot_rc.csv",
        discarded="results/df/unexpected_seqs/{sample}_unexpected.csv",
    conda:
        "../envs/jupyter.yaml"
    run:
        import pandas as pd

        unfiltered_df = pd.read_csv(input[0])
        unfiltered_df[unfiltered_df.Nham_aa <= config["random"]["Nham_aa_max"]].to_csv(
            output["filtered"], index=False
        )
        unfiltered_df[unfiltered_df.Nham_aa > config["random"]["Nham_aa_max"]].to_csv(
            output["discarded"], index=False
        )
