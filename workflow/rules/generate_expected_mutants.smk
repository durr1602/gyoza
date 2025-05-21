rule generate_mutants:
    input:
        layout=config["samples"]["path"],
        wtseqs=config["samples"]["wt"],
    output:
        config["samples"]["expected_mut"],
    message:
        f'Generating expected mutants based on the experimental design (codon mode = {config["codon"]["mode"]})'
    log:
        notebook="logs/notebooks/generate_mutants.ipynb",
    conda:
        "../envs/jupyter.yaml"
    notebook:
        "../notebooks/generate_mutants.py.ipynb"
