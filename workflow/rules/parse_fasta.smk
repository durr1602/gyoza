rule parse_fasta:
    input:
        fasta_files = expand(rules.vsearch.output, sample=samples),
        read_stats = rules.stats.output[0],
        expected_mutants = rules.generate_mutants.output[0]
    output:
        read_counts = 'results/df/readcounts.csv.gz',
        rc_filter_plot = report('results/graphs/rc_filter_plot.svg',
            '../report/rc_filter_plot.rst',
            category='1. Read filtering',
            labels={"figure": "1.1. Summary of filtered reads"}
        ),
        unexp_rc_plot = report('results/graphs/unexp_rc_plot.svg',
            '../report/unexp_rc_plot.rst',
            category='1. Read filtering',
            labels={"figure": "1.2. Read counts of unexpected variants"}
        )
    resources:
        mem_gb = 2, # > default to read csv.gz
        threads = 1,
        time = lambda _, attempt: f'00:{attempt}:00'
    message:
        "Parsing fasta files and comparing sequenced mutants with expectations..."
    log:
        notebook="logs/notebooks/parse_fasta.ipynb"
    conda:
        '../envs/jupyter_plotting.yaml'
    notebook:
        '../notebooks/parse_fasta.py.ipynb'
