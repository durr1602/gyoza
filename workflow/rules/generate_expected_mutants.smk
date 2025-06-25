rule generate_mutants:
    input:
        config["samples"]["wt"],
    output:
        f"{EXPMUT_PATH}/{{mutseq}}.csv.gz",
    message:
        f'Generating expected mutants based on the experimental design (codon mode = {config["codon"]["mode"]})'
    log:
        notebook="logs/0_generate_mutants/{mutseq}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/generate_mutants.py"
