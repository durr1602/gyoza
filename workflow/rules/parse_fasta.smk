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
        ),
        done = touch('results/done/parse_fasta.done')
    resources:
        mem_gb = lambda _, input, attempt: max(0.002*input.size_mb + (attempt-1)*0.002*input.size_mb, 1),
        threads = 1,
        time = lambda _, input, attempt: max((0.002*input.size_mb + (attempt-1)*0.002*input.size_mb).__ceil__(), 1)
    message:
        "Parsing fasta files and comparing sequenced mutants with expectations..."
    log:
        notebook="logs/notebooks/parse_fasta.ipynb"
    conda:
        '../envs/jupyter.yaml'
    notebook:
        '../notebooks/parse_fasta.py.ipynb'
