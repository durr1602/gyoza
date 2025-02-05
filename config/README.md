# Configuration of gyōza

## Project-specific files

A few files should be provided to properly analyze your data. What follows is the general procedure, however a toy dataset is provided in order to test the workflow. If you simply want to run the workflow with the toy dataset, enter the following: `snakemake --use-conda`.

### Sequencing data

Please provide the raw reads (forward and reverse) of your DMS sequencing data in the `config/reads` folder (or specify a different path in the main config file). The file names should be featured in the layout (see section below).

### Layout

Please provide a csv-formatted layout of your samples. The file should be named `layout.csv` and be located in the `config/project_files` folder. Here is [an example](project_files/layout.csv). The file should contain the following columns:
- Sample_name: the unique identifier for each of your samples. The sample name does not need to contain information about the timepoint or replicate, since these correspond to other columns
- R1: base name of the fastq file for forward (R1) reads (can be gzipped), including extension
- R2: base name of the fastq file for reverse (R2) reads (can be gzipped), including extension
- N_forward: the 5'-3' DNA sequence corresponding to the fixed region upstream of the mutated sequence
- N_reverse: the 5'-3' DNA sequence corresponding to the fixed region 5' of the mutated sequence on the reverse strand
- Mutated_seq: the unique identifier for the mutated DNA sequence, should be the same for all samples in which the same sequence was mutated
- Pos_start: starting position in the protein sequence. If you've mutated several regions/fragments in a coding gene, this position should refer to the full-length protein sequence
- Replicate: e.g. "R1"
- Timepoint: "T0", "T1", etc. Intermediate timepoints are optional.

Finally, additional columns can be added by the user to specify what makes this sample unique. These are referred to as "sample attributes" and could correspond to the genetic background, the fragment/region of the gene if it applies, the drug used for selection, etc. In summary, a "sample" is any unique combination of sample attributes + Replicate + Timepoint and should be associated to 2 fastq files, for the forward and reverse reads, respectively. Sample attributes = attributes related to Mutated_seq + optional attributes.

### WT DNA sequences

Please provide a csv-formatted list of WT DNA sequences. The file should be named `wt_seq.csv` and be located in the `config/project_files` folder. Here is [an example](project_files/wt_seq.csv). The file should contain **exactly** the two following columns:
- Mutated_seq: all possible values for the Mutated_seq flag from the layout
- WT_seq: corresponding WT DNA sequence

### Codon table

To prevent any typing mistake, the genetic code is imported from a [CoCoPUTs](https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs) table (which also features codon frequencies, although the workflow does not make use of this). [The one provided](project_files/ScerevisiaeTAXID559292_Cocoputs_codon_table.csv) corresponds to *Saccharomyces cerevisiae* TAXID 559292. Please edit the main [config file](config_file.yaml) if you ever need to specify a different genetic code. Any csv-formatted file with at least two columns ("codon" and "aminoacid") should do.

### Normalization with growth data

Normalization with growth data in the form of numbers of mitotic generations for each condition is **optional**. If you do want to normalize with growth data, please set the corresponding parameter to True (see section below on the main config file). In any case, a csv-formatted file will be **automatically generated** the first time the workflow is run (even if it is a dry run). Again, if normalization is set to True in the config, you will be prompted to edit the file to add the number of mitotic generations for each condition in the column 'Nb_gen'. Once the file is edited, re-run the workflow.

### Codon mode

Please specify the codon mode, meaning the type of degenerate codons you introduced at each position in the specified loci. Currently supported are: "NNN" (default value) or "NNK". Make sure you adapt the main [config file](config_file.yaml) if necessary.

## Main config file

The main config file is located [here](config_file.yaml). Please make sure to:
* select the samples to be processed (or leave 'all' if you want to process all samples)
* list your sample attributes
* replace all parameter values with the ones adapted for your project
* set the "normalize with growth" parameter to True if you want to normalize with growth
* set the "generate report" parameter to True if you want the HTML report to be automatically generated upon full completion of the workflow
* edit all directory/file paths if necessary

## Note on validation

Currently, all the following files are validated against a YAML schema to help spot formatting issues (misspelled column headers, missing mandatory properties, improper format, etc.): main config file, sample layout, file with WT DNA sequences, codon table, file with the number of mitotic generations.

## Technical configuration

The file containing technical config parameters to run the snakemake pipeline on HPC is [here](../profile/config.v8+.yaml). Apart from your email adress (please replace `<...>`), this file does not need to be modified too much, and flags added to the snakemake command line will supersede the default values specified in the file.
