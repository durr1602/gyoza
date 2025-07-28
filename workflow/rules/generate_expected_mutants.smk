def calc_mem(wildcards, input, attempt):
    # Manually calculate input filesize because input.size_mb gets it wrong in this specific case
    total_size = os.path.getsize(input[0])
    size_mb = total_size / 1024
    df = pd.read_csv(input[0], usecols=["codon_mode"])
    if (df.codon_mode.astype(str).str.count("x") == 1).any():
        factor = 1000
    else:
        factor = 1
    mem = max(0.05 * size_mb * factor * attempt, 2)
    return int(mem)


def calc_time(wildcards, input, attempt):
    # Manually calculate input filesize
    total_size = os.path.getsize(input[0])
    size_mb = total_size / 1024
    df = pd.read_csv(input[0], usecols=["codon_mode"])
    if (df.codon_mode.astype(str).str.count("x") == 1).any():
        factor = 500
    else:
        factor = 1
    alloc_time = max(0.08 * size_mb * factor * attempt, 2)
    return alloc_time


rule generate_mutants:
    input:
        config["samples"]["wt"],
    output:
        temp(f"{EXPMUT_PATH}/{{mutseq}}.csv.gz"),
    params:
        genetic_code=config["codon"]["table"],
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
