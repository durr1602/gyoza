def get_expmut_for_sample(wildcards):
    return str(EXPMUT_PATH / f"{sample_to_mutseq[wildcards.sample]}.csv.gz")


rule compare_to_sequencing:
    input:
        expected_mut_path=get_expmut_for_sample,
        readcount_path=rules.parse_fasta.output[0],
    output:
        observed="results/df/observed_seqs/{sample}_observed.csv",
        unexpected="results/df/unexpected_seqs/{sample}_unexpected.csv",
    params:
        readcount_level = RC_LEVEL,
    message:
        f"Comparing sequencing data to expected mutants"
    log:
        "logs/5_compare/compare_to_seq-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/compare_to_sequencing.py"
