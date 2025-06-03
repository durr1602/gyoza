rule process_read_counts:
    input:
        readcounts=lambda wildcards: expand(
            "results/df/annotated_readcounts/{sample}_annot_rc.csv",
            sample=grouped_samples_str[wildcards.group_key],
        ),
        nbgen=config["samples"]["generations"],
    output:
        selcoeffs="results/df/all_scores_{group_key}.csv",
        avg_scores="results/df/avg_scores_{group_key}.csv",
        hist_plot="results/graphs/hist_plot_{group_key}.svg",
        upset_plot="results/graphs/upset_plot_{group_key}.svg",
        timepoints_plot="results/graphs/timepoints_plot_{group_key}.svg",
        freq_df="results/df/distribution_freq/freq_{group_key}.csv",
        aa_df="results/df/agg_aa/aa_{group_key}.csv",
        #done=touch("results/done/process_read_counts.done"),
    message:
        "Processing read counts... converting to functional impact scores"
    log:
        notebook="logs/notebooks/process_read_counts_{group_key}.ipynb",
    conda:
        "../envs/jupyter.yaml"
    script:
        "../scripts/process_rc.py"
