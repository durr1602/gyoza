rule parse_fasta:
    input:
        fasta_files=expand(rules.vsearch_fastx_uniques.output, sample=samples),
        read_stats=rules.stats.output[0],
        expected_mutants=rules.annotate_mutants.output[0],
    output:
        read_counts="results/df/readcounts.csv.gz",
        rc_filter_plot=report(
            "results/graphs/rc_filter_plot.svg",
            "../report/rc_filter_plot.rst",
            category="1. Read filtering",
            labels={"figure": "1.1. Summary of filtered reads"},
        ),
        unexp_rc_plot=report(
            "results/graphs/unexp_rc_plot.svg",
            "../report/unexp_rc_plot.rst",
            category="1. Read filtering",
            labels={"figure": "1.2. Read counts of unexpected variants"},
        ),
        done=touch("results/done/parse_fasta.done"),
    message:
        "Parsing fasta files and comparing sequenced mutants with expectations..."
    log:
        notebook="logs/notebooks/parse_fasta.ipynb",
    conda:
        "../envs/jupyter.yaml"
    notebook:
        "../notebooks/parse_fasta.py.ipynb"
