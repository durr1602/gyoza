def extract_param_from_first_row(wildcards, input):
    df = pd.read_csv(input[0], nrows=1)
    if not df.loc[0, "WT"]:
        raise ValueError(
            f"Error.. First row of {input[0]} does not have WT == True. Cannot extract wild-type amino acid sequence."
        )

    return df.loc[0, "aa_seq"]


rule plot_read_counts:
    input:
        "results/df/annotated_readcounts/{sample}_annot_rc.csv",
    output:
        report(
            "results/graphs/heatmap_readcount_{sample}.svg",
            "../report/heatmap_rc.rst",
            category="1. Read filtering",
            subcategory="1.2. Heatmaps of raw read counts",
            labels={"figure": "{sample}"},
        ),
    params:
        wtaa=extract_param_from_first_row,
        level="codons",
    message:
        f"Plotting heatmaps of raw read counts.."
    log:
        "logs/7_plot/plot_rc-sample={sample}.log",
    conda:
        "../envs/main.yaml"
    script:
        "../scripts/plot_heatmaps.py"
