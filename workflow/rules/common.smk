##### Import libraries #####

import sys
import pandas as pd
from snakemake.utils import validate
from pathlib import Path
from collections import defaultdict
import warnings

##### Import and validate main config ####


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
print("Main config validated.")

##### Convert paths #####

READS_PATH = Path(config["reads"]["path"])
EXPMUT_PATH = Path(config["samples"]["expected_mut"])

##### Import and validate sample layout #####

layout_mandatory_cols = [
    "Sample_name",
    "R1",
    "R2",
    "N_forward",
    "N_reverse",
    "Mutated_seq",
    "Pos_start",
    "Replicate",
    "Timepoint",
]

layout_csv = pd.read_csv(config["samples"]["path"], dtype={"Replicate": str})

# Get timepoints other than T0
TIMEPOINTS = sorted(set(layout_csv["Timepoint"]) - {"T0"})

# Sanitize Replicate column
layout_csv["Replicate"] = layout_csv["Replicate"].fillna("").astype(str).str.strip()
validate(layout_csv, schema="../schemas/sample_layout.schema.yaml")

# Convert Report column to strict boolean
truthy = {"true", "1", "yes", "y", "on"}
layout_csv["Report"] = (
    layout_csv["Report"].fillna(False).astype(str).str.strip().str.lower().isin(truthy)
)

# Retrieve non-mandatory columns
layout_add_cols = [
    x for x in layout_csv.columns if x not in layout_mandatory_cols + ["Report"]
]

# Get sample / Mutated_seq mapping
sample_to_mutseq = dict(zip(layout_csv["Sample_name"], layout_csv["Mutated_seq"]))

sample_layout = layout_csv.set_index("Sample_name")
print("Sample layout validated.")

##### Validate sample attributes #####

if len(config["samples"]["attributes"]) > 0:
    for x in config["samples"]["attributes"]:
        if x not in layout_csv.columns:
            raise Exception(
                f"Error.. The sample attribute '{x}' does not feature as a column in the sample layout."
            )
        else:
            [
                warnings.warn(
                    f"The column '{x}' is not listed in your sample attributes"
                )
                for x in layout_add_cols
                if x not in config["samples"]["attributes"]
            ]
    TR_layout = layout_csv[["Sample_name"] + config["samples"]["attributes"]]
    grouped_samples = defaultdict(list)
    for _, row in TR_layout.iterrows():
        key = tuple(row[col] for col in config["samples"]["attributes"])
        grouped_samples[key].append(row["Sample_name"])
    print("Sample attributes imported.")
else:
    print("No sample attributes provided.")

##### Validate CSV file containing WT DNA sequences #####

if exists(config["samples"]["wt"]):
    wtseqs = pd.read_csv(config["samples"]["wt"])
    validate(wtseqs, schema="../schemas/wt_seqs.schema.yaml")
    mutseq_to_wtseq = dict(zip(wtseqs["Mutated_seq"], wtseqs["WT_seq"]))
    print("WT imported.")

##### Validate CSV file containing expected DNA sequences #####
# if exists(config["samples"]["expected_mut"]):
#    expmut = pd.read_csv(config["samples"]["expected_mut"])
#    validate(expmut, schema="../schemas/wt_seqs.schema.yaml")
#    print("Expected mutated sequences imported.")

##### Generate template CSV file to write the number of cellular generations between time points #####
# Note: At this time, this file is required to exist even if the user opts out of this normalization
# A template is generated with the column for the number of generations set to 1
# If the user opts out of normalization, this template will be used, dividing all scores by 1 (therefore no normalization)
# If the user opts in, a warning will notify the user that the template needs to be filled
# Once the column contains other values than 1 for every row, we'll use the data for normalization

if exists(config["samples"]["generations"]):
    nbgen = pd.read_csv(config["samples"]["generations"], dtype={"Replicate": str})
    validate(nbgen, schema="../schemas/nbgen.schema.yaml")
    if (config["normalize"]["with_gen"]) & ((nbgen.Nb_gen == 1).all(0)):
        raise Exception(
            f'>>Please fill in the file {config["samples"]["generations"]} with the number of cellular generations<<'
            + "\n>>(or deactivate this normalization in the main config file)<<"
        )
    elif config["normalize"]["with_gen"]:
        print("Ready to normalize with the provided numbers of cellular generations.")
    else:
        print("No normalization with cellular generations")
else:
    nbgen_temp = layout_csv[layout_csv.Timepoint != "T0"][
        config["samples"]["attributes"] + ["Replicate", "Timepoint"]
    ].drop_duplicates()
    nbgen_temp["Nb_gen"] = 1
    nbgen_temp.to_csv(config["samples"]["generations"], index=None)
    if config["normalize"]["with_gen"]:
        raise Exception(
            f'>>Please fill in the file {config["samples"]["generations"]} with the number of cellular generations<<'
            + "\n>>(or deactivate this normalization in the main config file)<<"
        )
    else:
        print("No normalization with cellular generations")

##### Validate codon table #####

codon_table = pd.read_csv(config["codon"]["table"], header=0)
validate(codon_table, schema="../schemas/codon_table.schema.yaml")
print("Codon table validated.")

##### Select samples to process #####

SAMPLES = sample_layout.sort_index().index

if config["samples"]["selection"] != "all":
    selection = [x for x in config["samples"]["selection"] if x in SAMPLES]
    if len(selection) == 0:
        raise Exception(
            "Error.. None of the samples you specified for processing feature in the sample layout."
        )
    elif len(selection) != len(config["samples"]["selection"]):
        statement = "Warning... at least one sample was misspelled when selecting samples to process in the config file\nWill continue with only the correctly spelled samples."
        warnings.warn(statement)
    else:
        SAMPLES = selection
        print("Selection of samples confirmed.")

MUTATED_SEQS = sorted(set(sample_to_mutseq[s] for s in SAMPLES))

grouped_samples = {
    key: [s for s in samples if s in SAMPLES]
    for key, samples in grouped_samples.items()
    if any(s in SAMPLES for s in samples)
}

##### Convert sample grouping wilcard <-> string #####


# Serialize tuple to string
def serialize_key(key):
    return "__".join(str(k) for k in key)


# Deserialize back to tuple
def deserialize_key(key_str):
    return tuple(key_str.split("__"))


# Map from string keys to sample lists
grouped_samples_str = {serialize_key(k): v for k, v in grouped_samples.items()}
GROUP_KEYS = list(grouped_samples_str.keys())

# Derive which groups should be reported
grouped_samples_for_report = [
    tuple(row[attr] for attr in config["samples"]["attributes"])
    for _, row in layout_csv[layout_csv["Report"]].iterrows()
]

# Ensure unique entries
seen = set()
grouped_samples_for_report = [
    x for x in grouped_samples_for_report if not (x in seen or seen.add(x))
]

# Filter out groups that don’t include any selected samples
grouped_samples_for_report = [
    group
    for group in grouped_samples_for_report
    if any(s in SAMPLES for s in grouped_samples.get(group, []))
]

# Serialize groups of samples to be included in the reported
REPORTED_GROUPS = [serialize_key(group) for group in grouped_samples_for_report]

# Derive individual samples to be reported
REPORTED_SAMPLES = sorted(
    {
        sample
        for group in grouped_samples_for_report
        for sample in grouped_samples.get(group, [])
        if sample in SAMPLES
    }
)


##### Cross groups x time points #####
GF = [(g, t) for g in GROUP_KEYS for t in TIMEPOINTS]
GF_REPORTED = [(g, t) for g in REPORTED_GROUPS for t in TIMEPOINTS]


##### Prepare HTML report #####
# Note: I've tried multiple approaches. Report cannot be reliably integrated
# in a dedicated rule (for DAG inclusion) because of the nested snakemake statements
# which unpredictably lead to filesystem errors.
# I went back to my initial approach (onsuccess hook) which felt hacky,
# but apparently is common practice (until something better comes)


def collect_graphs():
    graph_dir = Path("results/graphs")

    agg_graphs = [
        "rc_filter_plot.svg",
        "unexp_rc_plot.svg",
    ]

    group_specific_graphs = [f"heatmap_readcount_{s}.svg" for s in REPORTED_SAMPLES]

    if config["process_read_counts"]:
        agg_graphs += [
            "rc_var_plot.svg",
            "scoeff_violin_plot.svg",
            "replicates_heatmap_plot.svg",
            "replicates_plot.svg",
            "s_through_time_plot.svg",
        ]

        group_specific_graphs += (
            [f"hist_plot_{k}.svg" for k in REPORTED_GROUPS]
            + [f"upset_plot_{k}.svg" for k in REPORTED_GROUPS]
            + [f"timepoints_plot_{k}.svg" for k in REPORTED_GROUPS]
            + [f"heatmap_fitness_{k}_{t}.svg" for (k, t) in GF_REPORTED]
        )

    return [
        str(graph_dir / f)
        for f in agg_graphs + group_specific_graphs
        if (graph_dir / f).exists()
    ]


##### Specify final target #####


def get_target():
    targets = ["results/all_stats.csv"]
    targets += expand("results/graphs/heatmap_readcount_{sample}.svg", sample=SAMPLES)

    if config["process_read_counts"]:
        targets += expand("results/df/all_scores_{group_key}.csv", group_key=GROUP_KEYS)
        targets += expand(
            "results/graphs/heatmap_fitness_{group_key}_{t}.svg",
            group_key=GROUP_KEYS,
            t=TIMEPOINTS,
        )
        targets.append("results/graphs/rc_var_plot.svg")

    if config["qc"]["perform"]:
        targets.append("results/0_qc/multiqc.html")

    return targets
