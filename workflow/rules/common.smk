##### Setting rule agnostic env #####

conda:
    "../envs/global.yaml"

##### Import libraries #####

import pandas as pd
from snakemake.utils import validate
import warnings

##### Import main config ####

configfile: "config/config_file.yaml"

## Validate config
validate(config, schema="../schemas/config.schema.yaml")
print("Main config validated.")

##### Import and validate sample layout #####

layout_csv = pd.read_csv(config["samples"]["path"])
validate(layout_csv, schema="../schemas/sample_layout.schema.yaml")
sample_layout = layout_csv.set_index("Sample_name")
print("Sample layout validated.")

##### Validate TSV file containing WT DNA sequences #####
wtseqs = pd.read_csv(config["samples"]["wt"], sep='\t')
validate(wtseqs, schema="../schemas/wt_seqs.schema.yaml")
print("WT imported.")

##### Validate EXCEL file containing number of mitotic generations between time points #####
nbgen = pd.read_excel(config["samples"]["generations"])
validate(nbgen, schema="../schemas/nbgen.schema.yaml")
if (nbgen.Nb_gen == 1).all(0):
    warnings.warn("Numbers of mitotic generations all equal to 1 = No normalization with growth data.")
else:
    print("Ready to normalize with provided growth data.")

##### Validate sample attributes #####
for x in config["samples"]["attributes"]:
    if x not in layout_csv.columns:
        raise Exception(f"Error.. The sample attribute '{x}' does not feature as a column in the sample layout.")
    elif x not in nbgen.columns:
        raise Exception(f"Error.. The sample attribute '{x}' does not feature as a column in the EXCEL file with the number of mitotic generations.")
    else:
        pass
print('Sample attributes imported.')

##### Validate codon table #####
codon_table = pd.read_csv(config["codon"]["table"], header=0)
validate(codon_table, schema="../schemas/codon_table.schema.yaml")
print('Codon table validated.')

##### Select samples to process #####

samples = sample_layout.sort_index().index

if config["samples"]["selection"] != "all":
    selection = [x for x in config["samples"]["selection"] if x in samples]
    if len(selection) == 0:
        raise Exception("Error.. None of the samples you specified for processing feature in the sample layout.")
    elif len(selection) != len(config["samples"]["selection"]):
        statement = "Warning... at least one sample was misspelled when selecting samples to process in the config file\nWill continue with only the correctly spelled samples."
        warnings.warn(statement)
    else:
        samples = selection
        print('Selection of samples confirmed.')

##### Specify final target #####

def get_target():
    return ['results/df/selcoeffs.csv', 'results/0_qc/multiqc.html']
