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

##### Paths and other constant variables #####

PROJECT_DIR = Path(config["project"]["folder"])
READS_PATH = Path(config["reads"]["path"])
LAYOUT_PATH = PROJECT_DIR / "layout.csv"
GEN_CODE_PATH = PROJECT_DIR / "codon_table.csv"
WT_PATH = PROJECT_DIR / "wt_seq.csv"
EXPMUT_PATH = PROJECT_DIR / "expected_mut/"
NBGEN_PATH = PROJECT_DIR / "nbgen.csv"

SAMPLE_ATTR = config["project"]["sample_attributes"]
SCREEN_ATTR = config["project"]["screening_attributes"]

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

layout_csv = pd.read_csv(LAYOUT_PATH, dtype={"Replicate": str})

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
    if x not in SAMPLE_ATTR + SCREEN_ATTR:
        warnings.warn(
            f"Column {x} is not listed in your sample or screening attributes."
        )

if not SAMPLE_ATTR:
    raise ValueError(
        "Error.. Please specify at least one sample attribute (e.g. 'Mutated_seq')."
    )
else:
    for attr in SAMPLE_ATTR:
        if attr not in layout_csv.columns:
            raise Exception(f"Missing sample attribute column in the layout: {attr}.")

    print("Sample attributes imported.")

    # Initial sample grouping based on layout
    all_groups = defaultdict(list)

    # Map groups <-> samples by separating input (T0) and output
    # For each df, we need to convert the index to a tuple in case there's a single attribute
    tn = layout_csv.groupby(SAMPLE_ATTR + SCREEN_ATTR)["Sample_name"].agg(list)
    tn.index = tn.index.map(lambda x: (x,) if not isinstance(x, tuple) else x)
    t0 = (
        layout_csv[layout_csv.Timepoint == "T0"]
        .groupby(SAMPLE_ATTR)["Sample_name"]
        .agg(list)
    )
    t0.index = t0.index.map(lambda x: (x,) if not isinstance(x, tuple) else x)

    allg_df = tn.copy()

    for idx in allg_df.index:
        sample_key = idx[: len(SAMPLE_ATTR)]
        if sample_key in t0.index:
            allg_df[idx] = allg_df[idx] + t0[sample_key]
    all_groups = dict(allg_df)


##### Select samples to analyze/report #####


# Helper function to restrict to samples selected by the user
def select_samples(selection_column, sample_layout, all_groups):
    # Which sample_names are selected?
    selected_samples = set(sample_layout.index[sample_layout[selection_column]])

    groups = {}
    for group_key, samples in all_groups.items():
        # Keep only samples belonging to this group AND selected by the user
        selected_in_group = [s for s in samples if s in selected_samples]

        if selected_in_group:
            groups[group_key] = selected_in_group

    return groups


analyze_groups = select_samples("Analyze", sample_layout, all_groups)
report_groups = select_samples("Report", sample_layout, all_groups)

##### Merge all groups #####

if config["process_all_samples"]:
    report_groups = all_groups
    final_groups = all_groups
else:
    final_groups = {
        g: sorted(set(analyze_groups.get(g, []) + report_groups.get(g, [])))
        for g in set(analyze_groups) | set(report_groups)
    }

##### Final list of samples #####

SAMPLES = sorted({s for samples in final_groups.values() for s in samples})
REPORTED_SAMPLES = sorted({s for samples in report_groups.values() for s in samples})
MUTATED_SEQS = sorted(set(sample_to_mutseq[s] for s in SAMPLES))

T0_SAMPLES = [s for s in SAMPLES if sample_layout.loc[s, "Timepoint"] == "T0"]

if not T0_SAMPLES:
    raise WorkflowError(
        "Please select at least 1 T0 sample by writing Y in the Analyze column or in the Report column of the layout."
    )

print(f"{len(SAMPLES)} sample(s) selected for analysis.")


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

##### Get combinations of groups and output time points #####
# Keep only input/output pairs with at least one matching replicate

GT_WITH_OUTPUTS = []

for g, samples in final_groups.items():
    # Collect output samples
    outputs = [s for s in samples if sample_layout.loc[s, "Timepoint"] != "T0"]
    if not outputs:
        continue

    # Collect T0 replicates for matching
    t0_reps = {sample_layout.loc[s, "Replicate"] for s in T0_SAMPLES}

    # Collect if there's a matching T0 replicate
    outputs_with_matching_t0 = [
        s for s in outputs if sample_layout.loc[s, "Replicate"] in t0_reps
    ]

    if not outputs_with_matching_t0:
        continue

    # Get corresponding time point
    tp = sample_layout.loc[outputs_with_matching_t0[0], "Timepoint"]

    GT_WITH_OUTPUTS.append((serialize_key(g), tp))

if not GT_WITH_OUTPUTS:
    raise WorkflowError(
        "Error.. Please select at least one pair of matching input/output replicates."
    )

ATTR_GROUPS_WITH_OUTPUTS = sorted({g for (g, tp) in GT_WITH_OUTPUTS})
REPORTED_GROUPS_WITH_OUTPUTS = sorted(
    {g for (g, tp) in GT_WITH_OUTPUTS if g in REPORTED_GROUPS}
)

##### Validate CSV file containing WT DNA sequences #####
# Required only for 'codon' and 'random' designs
# If avaible, we retrieve the WT from there to be able to annotate mutants

mutseq_to_wtseq = {}

if exists(WT_PATH):
    wtseqs = pd.read_csv(WT_PATH)
    validate(wtseqs, schema="../schemas/wt_seqs.schema.yaml")

    if set(wtseqs.Mutated_seq.unique()).isdisjoint(set(MUTATED_SEQS)):
        raise WorkflowError(
            f"Error.. None of the Mutated_seq values in {WT_PATH} match those in {LAYOUT_PATH}"
        )

    wtseqs["WT_seq"] = wtseqs["WT_seq"].str.upper()
    mutseq_to_wtseq = dict(zip(wtseqs["Mutated_seq"], wtseqs["WT_seq"]))
    print("WT imported.")

##### Validate CSV files containing expected DNA sequences #####
# Note: WT CSV is not required for 'provided' design
# For 'provided' and 'random' designs, we get the WT from the list of expected mutants

expmut_mutseqs = []

for f in EXPMUT_PATH.glob("*.csv.gz"):
    expmut = pd.read_csv(f)
    validate(expmut, schema="../schemas/exp_mut.schema.yaml")
    if expmut["Mutated_seq"].nunique() != 1:
        raise ValueError(f"Error.. Multiple 'Mutated_seq' values in {f.name}")
    mutseq = expmut.at[0, "Mutated_seq"]
    expmut_mutseqs.append(mutseq)
    if expmut["WT_seq"].nunique() != 1:
        raise ValueError(f"Error.. Multiple 'WT_seq' values in {f.name}")
    wtseq = expmut.at[0, "WT_seq"].upper()
    mutseq_to_wtseq[mutseq] = wtseq
    if not (expmut["nt_seq"].astype(str).str.len() == len(wtseq)).all():
        raise ValueError(
            f"Not all 'nt_seq's have the same length as 'WT_seq' in {f.name}"
        )
    print(f"Imported expectant mutants of {mutseq}.")
else:
    expmut_mutseqs = MUTATED_SEQS

if set(expmut_mutseqs).isdisjoint(set(MUTATED_SEQS)):
    raise WorkflowError(
        f"Error.. None of the Mutated_seq values imported from files in {EXPMUT_PATH} match those in {LAYOUT_PATH}"
    )

##### Validate codon table #####

codon_table = pd.read_csv(GEN_CODE_PATH, header=0)
validate(codon_table, schema="../schemas/codon_table.schema.yaml")
print("Codon table validated.")
codon_table["aminoacid"] = codon_table["aminoacid"].str.upper()
codon_table["codon"] = codon_table["codon"].str.upper()
GEN_CODE = dict(zip(codon_table["codon"], codon_table["aminoacid"]))


# Define function to translate any DNA sequence
def get_aa_seq(nt, codon_dict):
    r"""Translates nucleotide sequence to amino acid sequence from codon dict.

    Parameters
    ----------
    nt : str
        DNA sequence (length should be a multiple of 3).
    codon_dict : dict
        Codon table associating codons to amino acid residues.

    Returns
    -------
    str

    Raises
    ------
    ValueError
        If the length of `nt` is not a multiple of 3.
    """
    if len(nt) % 3 != 0:
        raise ValueError(
            f"Error.. the length of the DNA sequence is not a multiple of 3."
        )

    nt_codons = [nt[i : i + 3] for i in range(0, len(nt), 3)]
    aa = "".join([codon_dict.get(x) for x in nt_codons])

    return aa


# Map WT amino acid sequence for each mutated locus
mutseq_to_wtaa = {
    mutseq: get_aa_seq(wtseq, GEN_CODE) for mutseq, wtseq in mutseq_to_wtseq.items()
}

# Map WT amino acid sequence for each group
group_to_wtaa = {
    group_key: get_aa_seq(mutseq_to_wtseq[sample_to_mutseq[samples[0]]], GEN_CODE)
    for group_key, samples in final_groups_str.items()
}

##### Generate template CSV file to write the number of cellular generations between time points #####
# Note: At this time, this file is required to exist even if the user opts out of this normalization
# A template is generated with the column for the number of generations set to 1
# If the user opts out of normalization, this template will be used, dividing all scores by 1 (therefore no normalization)
# If the user opts in, a warning will notify the user that the template needs to be filled
# Once the column contains other values than 1 for every row, we'll use the data for normalization

required_rows = layout_csv[
    layout_csv["Sample_name"].isin(SAMPLES) & (layout_csv["Timepoint"] != "T0")
][SAMPLE_ATTR + SCREEN_ATTR + ["Replicate", "Timepoint"]].drop_duplicates()

if exists(NBGEN_PATH):
    nbgen = pd.read_csv(NBGEN_PATH, dtype={"Replicate": str})
    validate(nbgen, schema="../schemas/nbgen.schema.yaml")
    if (config["normalize_with_gen"]) & ((nbgen.Nb_gen == 1).any()):
        raise Exception(
            f">>Please fill in the file {NBGEN_PATH} with the number of cellular generations<<\n"
            ">>(or deactivate this normalization in the main config file)<<"
        )
    elif config["normalize_with_gen"]:
        # Find missing rows from existing file
        merged = required_rows.merge(
            nbgen,
            on=SAMPLE_ATTR + SCREEN_ATTR + ["Replicate", "Timepoint"],
            how="left",
            indicator=True,
        )
        missing_rows = merged[merged["_merge"] == "left_only"].drop(columns=["_merge"])
        missing_rows["Nb_gen"] = 1

        # Append to existing file
        if not missing_rows.empty:
            print(f"Adding {len(missing_rows)} missing row(s) to {NBGEN_PATH}")
            nbgen = pd.concat([nbgen, missing_rows], ignore_index=True)
            nbgen.to_csv(NBGEN_PATH, index=False)

        # Additional check to make sure all rows from selection are filled properly
        # This is done to prevent bothering the user with filling data for non currently selected samples
        selected_rows = nbgen.merge(
            required_rows,
            on=SAMPLE_ATTR + SCREEN_ATTR + ["Replicate", "Timepoint"],
            how="inner",
        )

        if (selected_rows["Nb_gen"] == 1).any():
            raise Exception(
                f">> Please fill in the file {NBGEN_PATH} with the number of cellular generations <<\n"
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
    nbgen_temp.to_csv(NBGEN_PATH, index=None)
    if config["normalize_with_gen"]:
        raise Exception(
            f">> Please fill in {NBGEN_PATH} with the number of cellular generations <<\n"
            ">> Or disable this normalization in the config <<"
        )
    else:
        print("No normalization with cellular generations")


##### Helper functions for dynamic allocation of resources #####
def calc_mem(wildcards, input, attempt):
    # Manually calculate input filesize because input.size_mb gets it wrong in this specific case
    total_size = os.path.getsize(input[0])
    size_mb = total_size / 1024
    df = pd.read_csv(input[0], usecols=["codon_mode"])
    if (df.codon_mode.astype(str).str.count("x") == 1).any():
        factor = 1000
    else:
        factor = 1
    mem = max(0.05 * size_mb * factor * attempt, 2)
    return int(mem)


def calc_time(wildcards, input, attempt):
    # Manually calculate input filesize
    total_size = os.path.getsize(input[0])
    size_mb = total_size / 1024
    df = pd.read_csv(input[0], usecols=["codon_mode"])
    if (df.codon_mode.astype(str).str.count("x") == 1).any():
        factor = 500
    else:
        factor = 1
    alloc_time = max(0.08 * size_mb * factor * attempt, 2)
    return alloc_time


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
            + [f"heatmap_fitness_{k}_{t}.svg" for (k, t) in GT_WITH_OUTPUTS]
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
        if config["perform_qc"] and qc_path.exists():
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
    targets = ["results/df/all_stats.csv"]
    targets += expand("results/graphs/heatmap_readcount_{sample}.svg", sample=SAMPLES)

    if config["process_read_counts"]:
        targets.append(["results/df/all_scores.csv", "results/graphs/rc_var_plot.svg"])
        targets += expand(
            "results/graphs/heatmap_fitness_{group_key}_{t}.svg",
            zip,
            group_key=[g for g, t in GT_WITH_OUTPUTS],
            t=[t for g, t in GT_WITH_OUTPUTS],
        )

    if config["perform_qc"]:
        targets.append("results/0_qc/multiqc.html")

    return targets
