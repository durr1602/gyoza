def calc_mem(wildcards, input, attempt):
    # Manually calculate input filesize because input.size_mb gets it wrong in this specific case
    total_size = 0
    for f in input:
        total_size += os.path.getsize(f)
    size_mb = total_size / 1024
    factor = 1000 if config["codon"]["mode"].count("x") == 1 else 1
    mem = max(0.05 * size_mb * factor * attempt, 2)
    return int(mem)


def calc_time(wildcards, input, attempt):
    # Manually calculate input filesize
    total_size = 0
    for f in input:
        total_size += os.path.getsize(f)
    size_mb = total_size / 1024
    factor = 500 if config["codon"]["mode"].count("x") == 1 else 1
    alloc_time = max(0.08 * size_mb * factor * attempt, 2)
    return alloc_time


rule generate_mutants:
    input:
        config["samples"]["wt"],
    output:
        temp(f"{EXPMUT_PATH}/{{mutseq}}.csv.gz"),
    resources:
        mem_gb=calc_mem,
        time=calc_time,
    message:
        f'Generating expected mutants based on the experimental design (codon mode = {config["codon"]["mode"]})'
    log:
        notebook="logs/0_generate_mutants/{mutseq}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/generate_mutants.py"
