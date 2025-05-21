rule annotate_mutants:
    input:
        layout=config["samples"]["path"],
        exp_mut=config["samples"]["expected_mut"],
    output:
        "results/df/master_layout.csv.gz",
    message:
        f"Annotating expected mutants (number of mutations, alternative codons, etc)"
    log:
        notebook="logs/notebooks/annotate_mutants.ipynb",
    conda:
        "../envs/jupyter.yaml"
    notebook:
        "../notebooks/annotate_mutants.py.ipynb"
