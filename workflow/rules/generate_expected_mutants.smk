rule generate_mutants:
    input:
        WT_PATH,
    output:
        temp(f"{EXPMUT_PATH}/{{mutseq}}.csv.gz"),
    params:
        genetic_code=GEN_CODE_PATH,
    resources:
        mem_gb=calc_mem,
        time=calc_time,
    message:
        f"Generating expected mutants.."
    log:
        notebook="logs/0_generate_mutants/{mutseq}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/generate_mutants.py"
