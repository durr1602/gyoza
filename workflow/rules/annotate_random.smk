rule annotate_random:
    input:
        rules.parse_fasta.output
    output:
        annot_rc = "results/df/annotated_readcounts/{sample}_annot_rc.csv",
        indels = "results/df/indels/{sample}_indels.csv",
    message:
        f"Annotating mutants observed in sequencing data"
    log:
        "logs/6_annotate/annotate_random-sample={sample}.log",
    conda:
        "../envs/jupyter.yaml"
    script:
        "../scripts/annotate_mutants.py"

rule write_empty_unexpected:
    output:
        "results/df/unexpected_seqs/{sample}_unexpected.csv"
    conda:
        "../envs/jupyter.yaml"
    run:
        import pandas as pd
        pd.DataFrame(columns=["Sample_name","Mutated_seq","WT_seq","nt_seq","readcount"]).to_csv(output[0], index=False)