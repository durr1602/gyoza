##### Import libraries #####

import sys
import pandas as pd
from snakemake.utils import validate
from pathlib import Path
from collections import defaultdict
import warnings

##### Import and validate main config #####


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
    "Analyze",
    "Report",
]

layout_csv = pd.read_csv(config["samples"]["path"], dtype={"Replicate": str})

# Sanitize Replicate column
layout_csv["Replicate"] = layout_csv["Replicate"].fillna("").astype(str).str.strip()
validate(layout_csv, schema="../schemas/sample_layout.schema.yaml")

# Convert Report column to strict boolean
truthy = {"true", "t", "yes", "y", "ok", "1"}
for col in ["Analyze", "Report"]:
    layout_csv[col] = (
        layout_csv[col].fillna("").astype(str).str.strip().str.lower().isin(truthy)
    )

# Retrieve non-mandatory columns
layout_add_cols = [x for x in layout_csv.columns if x not in layout_mandatory_cols]

# Get sample / Mutated_seq mapping
sample_to_mutseq = dict(zip(layout_csv["Sample_name"], layout_csv["Mutated_seq"]))

sample_layout = layout_csv.set_index("Sample_name").sort_index()
print("Sample layout validated.")

##### Validate sample attributes and group samples #####

for x in layout_add_cols:
    if x not in config["samples"]["attributes"]:
        warnings.warn(f"Column {x} is not listed in your sample attributes.")

if not config["samples"]["attributes"]:
    raise ValueError("Error.. Please specify at least one sample attribute.")
else:
    for attr in config["samples"]["attributes"]:
        if attr not in layout_csv.columns:
            raise Exception(f"Missing sample attribute column in the layout: {attr}.")

    TR_layout = layout_csv[["Sample_name"] + config["samples"]["attributes"]]
    # Initial sample grouping based on layout
    all_groups = defaultdict(list)
    for _, row in TR_layout.iterrows():
        key = tuple(row[col] for col in config["samples"]["attributes"])
        all_groups[key].append(row["Sample_name"])
    print("Sample attributes imported.")


##### Select samples to analyze/report #####


# Helper function to rescue T0 + matching output timepoint replicates
def select_samples(selection_column, sample_layout, all_groups):
    selected_samples = sample_layout[sample_layout[selection_column]].index.tolist()
    groups = {}
    for group, samples in all_groups.items():
        selected_in_group = [s for s in samples if s in selected_samples]
        if not selected_in_group:
            continue
        selected_timepoints = {
            sample_layout.loc[s, "Timepoint"]
            for s in selected_in_group
            if sample_layout.loc[s, "Timepoint"] != "T0"
        }
        groups[group] = sorted(
            [
                s
                for s in samples
                if sample_layout.loc[s, "Timepoint"] == "T0"
                or sample_layout.loc[s, "Timepoint"] in selected_timepoints
            ]
        )
    return groups


analyze_groups = select_samples("Analyze", sample_layout, all_groups)
report_groups = select_samples("Report", sample_layout, all_groups)

##### Merge all groups #####

final_groups = defaultdict(list)
for group, samples in analyze_groups.items():
    final_groups[group].extend(samples)
for group, samples in report_groups.items():
    final_groups[group].extend(samples)
final_groups = {group: sorted(set(samples)) for group, samples in final_groups.items()}


##### Convert sample grouping wilcard <-> string #####


# Serialize tuple to string
def serialize_key(key):
    return "__".join(str(k) for k in key)


# Deserialize back to tuple
def deserialize_key(key_str):
    return tuple(key_str.split("__"))


# Serialize sample groups
final_groups_str = {
    serialize_key(group): samples for group, samples in final_groups.items()
}
ATTR_GROUPS = list(final_groups_str.keys())

REPORTED_GROUPS = [serialize_key(g) for g in report_groups if g in final_groups]

# Map back position position offset for each group
pos_offset_by_group = {
    group_key: sample_layout.loc[samples[0], "Pos_start"]
    for group_key, samples in final_groups_str.items()
}

##### Determine groups with output timepoints #####

groups_with_output_timepoints = {}
for group, samples in final_groups.items():
    timepoints = {sample_layout.loc[s, "Timepoint"] for s in samples}
    if "T0" in timepoints and any(tp != "T0" for tp in timepoints):
        groups_with_output_timepoints[group] = samples

ATTR_GROUPS_WITH_OUTPUTS = [
    serialize_key(group) for group in groups_with_output_timepoints
]
REPORTED_GROUPS_WITH_OUTPUTS = [
    g for g in ATTR_GROUPS_WITH_OUTPUTS if g in REPORTED_GROUPS
]

##### Get combinations of groups and output time points #####

GT_WITH_OUTPUTS = []
for group, samples in groups_with_output_timepoints.items():
    gkey = serialize_key(group)
    timepoints = sorted(
        {
            sample_layout.loc[s, "Timepoint"]
            for s in samples
            if sample_layout.loc[s, "Timepoint"] != "T0"
        }
    )
    GT_WITH_OUTPUTS.extend([(gkey, t) for t in timepoints])

GT_REPORTED = [(g, t) for (g, t) in GT_WITH_OUTPUTS if g in REPORTED_GROUPS]

##### Final list of samples #####

SAMPLES = sorted({s for samples in final_groups.values() for s in samples})
REPORTED_SAMPLES = sorted({s for samples in report_groups.values() for s in samples})
MUTATED_SEQS = sorted(set(sample_to_mutseq[s] for s in SAMPLES))

if not SAMPLES:
    raise Exception("No samples marked for analysis in the layout.")

print(f"{len(SAMPLES)} sample(s) selected for analysis.")

##### Validate CSV file containing WT DNA sequences #####
# Required only for 'codon' and 'random' designs
# If avaible, we retrieve the WT from there to be able to annotate mutants

mutseq_to_wtseq = {}

if exists(config["samples"]["wt"]):
    wtseqs = pd.read_csv(config["samples"]["wt"])
    validate(wtseqs, schema="../schemas/wt_seqs.schema.yaml")
    mutseq_to_wtseq = dict(zip(wtseqs["Mutated_seq"], wtseqs["WT_seq"]))
    print("WT imported.")

##### Validate CSV files containing expected DNA sequences #####
# Note: WT CSV is not required for 'provided' design
# For 'provided' and 'random' designs, we get the WT from the list of expected mutants

for f in EXPMUT_PATH.glob("*.csv.gz"):
    expmut = pd.read_csv(f)
    validate(expmut, schema="../schemas/wt_seqs.schema.yaml")
    mutseq = expmut.at[0, "Mutated_seq"]
    wtseq = expmut.at[0, "WT_seq"].upper()
    mutseq_to_wtseq[mutseq] = wtseq
    print(f"Imported expectant mutants of {mutseq}.")

##### Validate codon table #####

codon_table = pd.read_csv(config["codon"]["table"], header=0)
validate(codon_table, schema="../schemas/codon_table.schema.yaml")
print("Codon table validated.")

##### Generate template CSV file to write the number of cellular generations between time points #####
# Note: At this time, this file is required to exist even if the user opts out of this normalization
# A template is generated with the column for the number of generations set to 1
# If the user opts out of normalization, this template will be used, dividing all scores by 1 (therefore no normalization)
# If the user opts in, a warning will notify the user that the template needs to be filled
# Once the column contains other values than 1 for every row, we'll use the data for normalization

required_rows = layout_csv[
    layout_csv["Sample_name"].isin(SAMPLES) & (layout_csv["Timepoint"] != "T0")
][config["samples"]["attributes"] + ["Replicate", "Timepoint"]].drop_duplicates()

if exists(config["samples"]["generations"]):
    nbgen = pd.read_csv(config["samples"]["generations"], dtype={"Replicate": str})
    validate(nbgen, schema="../schemas/nbgen.schema.yaml")
    if (config["normalize"]["with_gen"]) & ((nbgen.Nb_gen == 1).any()):
        raise Exception(
            f">>Please fill in the file {config['samples']['generations']} with the number of cellular generations<<\n"
            ">>(or deactivate this normalization in the main config file)<<"
        )
    elif config["normalize"]["with_gen"]:
        # Find missing rows from existing file
        merged = required_rows.merge(
            nbgen,
            on=config["samples"]["attributes"] + ["Replicate", "Timepoint"],
            how="left",
            indicator=True,
        )
        missing_rows = merged[merged["_merge"] == "left_only"].drop(columns=["_merge"])
        missing_rows["Nb_gen"] = 1

        # Append to existing file
        if not missing_rows.empty:
            print(
                f"Adding {len(missing_rows)} missing row(s) to {config['samples']['generations']}"
            )
            nbgen = pd.concat([nbgen, missing_rows], ignore_index=True)
            nbgen.to_csv(config["samples"]["generations"], index=False)

        # Additional check to make sure all rows from selection are filled properly
        # This is done to prevent bothering the user with filling data for non currently selected samples
        selected_rows = nbgen.merge(
            required_rows,
            on=config["samples"]["attributes"] + ["Replicate", "Timepoint"],
            how="inner",
        )

        if (selected_rows["Nb_gen"] == 1).any():
            raise Exception(
                f">> Please fill in the file {config['samples']['generations']} with the number of cellular generations <<\n"
                ">>(or deactivate this normalization in the main config file)<<"
            )
        else:
            print(
                "Ready to normalize with the provided numbers of cellular generations."
            )
    else:
        print("No normalization with cellular generations")
else:
    nbgen_temp = required_rows
    nbgen_temp["Nb_gen"] = 1
    nbgen_temp.to_csv(config["samples"]["generations"], index=None)
    if config["normalize"]["with_gen"]:
        raise Exception(
            f">> Please fill in {config['samples']['generations']} with the number of cellular generations <<\n"
            ">> Or disable this normalization in the config <<"
        )
    else:
        print("No normalization with cellular generations")


##### Prepare HTML report #####
# Note: I've tried multiple approaches. Report cannot be reliably integrated
# in a dedicated rule (for DAG inclusion) because of the nested snakemake statements
# which unpredictably lead to filesystem errors.
# I went back to my initial approach (onsuccess hook) which felt hacky,
# but apparently is common practice (until something better comes)


def collect_graphs():
    graph_dir = Path("results/graphs")

    agg_graphs = ["rc_filter_plot.svg", "unexp_rc_plot.svg"]
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
            + [f"heatmap_fitness_{k}_{t}.svg" for (k, t) in GT_REPORTED]
        )

    return [
        str(graph_dir / f)
        for f in agg_graphs + group_specific_graphs
        if (graph_dir / f).exists()
    ]


def generate_report():
    try:
        graphs = collect_graphs()

        # Add QC report if present
        qc_path = Path("results/0_qc/multiqc.html")
        if config["qc"]["perform"] and qc_path.exists():
            graphs.append(str(qc_path))

        if graphs:
            # Inject config (if specified on the command line)
            config_arg = None
            args = sys.argv
            if "--configfile" in args:
                idx = args.index("--configfile")
                if idx + 1 < len(args):
                    config_arg = args[idx + 1]

            configfile_str = f"--configfile {config_arg}" if config_arg else ""

            # Switch to .zip if too many graphs collected
            ext = "html" if len(graphs) <= 30 else "zip"

            # CSS Style sheet
            css = "report-stylesheet"

            report_cmd = f"snakemake {' '.join(str(f) for f in graphs)} --report results/report.{ext} --{css} config/style/{css}.css {configfile_str}"
            shell(report_cmd)
        else:
            print(">> No graphs found. Skipping report generation.")

    except Exception as e:
        print(f"Report generation failed: {e}")


##### Workflow targets #####


def get_target():
    targets = ["results/all_stats.csv"]
    targets += expand("results/graphs/heatmap_readcount_{sample}.svg", sample=SAMPLES)

    if config["process_read_counts"]:
        targets += expand(
            "results/df/all_scores_{group_key}.csv",
            group_key=ATTR_GROUPS_WITH_OUTPUTS,
        )
        targets += expand(
            "results/graphs/heatmap_fitness_{group_key}_{t}.svg",
            zip,
            group_key=[g for g, t in GT_WITH_OUTPUTS],
            t=[t for g, t in GT_WITH_OUTPUTS],
        )
        targets.append("results/graphs/rc_var_plot.svg")

    if config["qc"]["perform"]:
        targets.append("results/0_qc/multiqc.html")

    return targets
