rule generate_mutants:
    params:
        layout = config["samples"]["path"],
        seqs = config["samples"]["wt"], # Projet-specific file containing the wild-type sequences
        codon_table = config["codon"]["table"], # Projet-specific file containing the genetic code
        codon_mode = config["codon"]["mode"] # Project-specific parameter to specify which degenerate codons were introduced
    output:
        'results/df/master_layout.csv.gz'
    resources:
        threads = 1,
        time = lambda _, attempt: f'00:{attempt}:00'
    message:
        "Generating expected mutants based on the experimental design (codon mode = {params.codon_mode})"
    log:
        notebook="logs/notebooks/generate_mutants.ipynb"
    conda:
        '../envs/jupyter_basic.yaml'
    notebook:
        '../notebooks/generate_mutants.py.ipynb'
