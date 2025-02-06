rule generate_mutants:
    output:
        'results/df/master_layout.csv.gz'
    resources:
        threads = 1,
        time = lambda _, attempt: f'00:{attempt}:00'
    message:
        f'Generating expected mutants based on the experimental design (codon mode = {config["codon"]["mode"]})'
    log:
        notebook="logs/notebooks/generate_mutants.ipynb"
    conda:
        '../envs/jupyter_basic.yaml'
    notebook:
        '../notebooks/generate_mutants.py.ipynb'
