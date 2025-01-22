rule process_read_counts:
    input:
        rules.parse_fasta.output.read_counts
    params:
        sample_attributes = config['samples']['attributes'],
        exp_rc_per_var = config['rc_aims']['exp_rc_per_var'],
        rc_threshold = config['filter']['rc_threshold'],
        nb_gen = config['samples']['generations']
    output:
        selcoeffs = 'results/df/selcoeffs.csv',
        hist_plot = report('results/graphs/hist_plot.svg',
            '../report/hist_plot.rst',
            category='2. Read processing',
            labels={"figure": "2.1. Raw read count per variant"}
        ),
        upset_plot = report('results/graphs/upset_plot.svg',
            '../report/upset_plot.rst',
            category='2. Read processing',
            labels={"figure": "2.2. Overlap across time points and replicates"}
        ),
        rc_var_plot = report('results/graphs/rc_var_plot.svg',
            '../report/rc_var_plot.rst',
            category="2. Read processing",
            labels={"figure": "2.3. Distribution of allele frequencies"}
        ),
        timepoints_plot = report('results/graphs/timepoints_plot.svg',
            '../report/timepoints_plot.rst',
            category="3. Functional impact",
            labels={"figure": "3.3. Correlation between time points"}
        ),
        scoeff_violin_plot = report('results/graphs/scoeff_violin_plot.svg',
            '../report/scoeff_violin_plot.rst',
            category="3. Functional impact",
            labels={"figure": "3.1. Distribution of functional impact scores"}
        ),
        s_through_time_plot = report('results/graphs/s_through_time_plot.svg',
            '../report/s_through_time_plot.rst',
            category="3. Functional impact",
            labels={"figure": "3.4. Functional impact over time"}
        ),
        replicates_plot = report('results/graphs/replicates_plot.svg',
            '../report/replicates_plot.rst',
            category="3. Functional impact",
            labels={"figure": "3.2. Correlation between replicates"}
        )
    resources:
        mem_gb = 2, # > default to read csv.gz
        threads = 1,
        time = "00:02:00"
    message:
        "Processing read counts... converting to functional impact scores"
    log:
        notebook="logs/notebooks/process_read_counts.ipynb"
    conda:
        '../envs/jupyter_plotting.yaml'
    notebook:
        '../notebooks/process_read_counts.py.ipynb'
