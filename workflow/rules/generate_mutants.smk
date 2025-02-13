rule generate_mutants:
    output:
        'results/df/master_layout.csv.gz'
    resources:
        mem_gb = lambda _, totalNbCodons, attempt: max(0.005*totalNbCodons + (attempt-1)*0.005*totalNbCodons, 1),
        threads = 1,
        time = lambda _, totalNbCodons, attempt: max(0.01*totalNbCodons + (attempt-1)*0.01*totalNbCodons, 1)
    message:
        f'Generating expected mutants based on the experimental design (codon mode = {config["codon"]["mode"]})'
    log:
        notebook="logs/notebooks/generate_mutants.ipynb"
    conda:
        '../envs/jupyter.yaml'
    notebook:
        '../notebooks/generate_mutants.py.ipynb'
