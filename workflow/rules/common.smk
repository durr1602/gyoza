##### Import libraries #####

import pandas as pd
from snakemake.utils import validate
from pathlib import Path
import warnings

##### Import and validate main config ####


configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")
print("Main config validated.")

##### Convert paths #####

READS_PATH = Path(config['reads']['path'])
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
layout_csv = pd.read_csv(config["samples"]["path"])
validate(layout_csv, schema="../schemas/sample_layout.schema.yaml")
layout_add_cols = [x for x in layout_csv.columns if x not in layout_mandatory_cols]
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
#if exists(config["samples"]["expected_mut"]):
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
    nbgen = pd.read_csv(config["samples"]["generations"])
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

##### Specify final target #####

def get_target():
    targets = ["results/df/all_scores.csv"]
    if config["qc"]["perform"]:
        targets.append("results/0_qc/multiqc.html")
    return targets
