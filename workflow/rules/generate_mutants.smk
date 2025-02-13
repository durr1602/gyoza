rule generate_mutants:
    input:
        layout = config["samples"]["path"],
        wtseqs = config["samples"]["wt"]
    output:
        'results/df/master_layout.csv.gz'
    resources:
        mem_gb = lambda _, input, attempt: max(50*input.size_mb + (attempt-1)*50*input.size_mb, 2),
        threads = 1,
        time = lambda _, input, attempt: max(80*input.size_mb + (attempt-1)*80*input.size_mb, 1)
    message:
        f'Generating expected mutants based on the experimental design (codon mode = {config["codon"]["mode"]})'
    log:
        notebook="logs/notebooks/generate_mutants.ipynb"
    conda:
        '../envs/jupyter.yaml'
    notebook:
        '../notebooks/generate_mutants.py.ipynb'
