# Configuration of gyōza

## Project-specific files

A few files should be provided to properly analyze your data. What follows is the general procedure, however a toy dataset is provided in order to test the workflow. If you simply want to run the workflow with the toy dataset, enter the following: `snakemake`.

Unless specified otherwise, most file paths and parameters can be modified in the main [config file](config.yaml).

### Sequencing data

Please provide the raw reads (forward and reverse) of your DMS sequencing data in the `config/reads` folder (or specify a different path in the main config file). The file names should be featured in the layout (see section below).

### Layout

Please provide a csv-formatted layout of your samples. The file should be named `layout.csv` and be located in the `config/project_files` folder. Here is [an example](project_files/layout.csv). The file should contain the following columns:
- Sample_name: the unique identifier for each of your samples. The sample name does not need to contain information about the timepoint or replicate, since these correspond to other columns
- R1: base name of the fastq file for forward (R1) reads (can be gzipped), including extension
- R2: base name of the fastq file for reverse (R2) reads (can be gzipped), including extension
- N_forward: the 5'-3' DNA sequence corresponding to the fixed region upstream of the mutated sequence or anything that can be used as -g flag with cutadapt (including complex patterns such as 'NNATG;optional...ATG', in which case do not forget the single quotes)
- N_reverse: the 5'-3' DNA sequence corresponding to the fixed region 5' of the mutated sequence on the reverse strand or anything that can be used as -G flag with cutadapt (same requirements as above)
- Mutated_seq: the unique identifier for the mutated DNA sequence, should be the same for all samples in which the same sequence was mutated
- Pos_start: starting position in the protein sequence. If you've mutated several regions/fragments in a coding gene, this position should refer to the full-length protein sequence
- Replicate: e.g. "R1"
- Timepoint: "T0", "T1", etc. Intermediate timepoints are optional.
- Report: "yes" (or a different truthy value) to include the sample in the HTML report.

Finally, additional columns can be added by the user to specify what makes this sample unique. These are referred to as "sample attributes" and could correspond to the genetic background, the fragment/region of the gene if it applies, the drug used for selection, etc. In summary, a "sample" is any unique combination of sample attributes + Replicate + Timepoint and should be associated to 2 fastq files, for the forward and reverse reads, respectively. Sample attributes = attributes related to Mutated_seq + optional attributes.

### WT DNA sequences

If you want gyōza to automatically generate all expected sequences based on a codon mode, please provide a csv-formatted list of WT DNA sequences. The file should be named `wt_seq.csv` and be located in the `config/project_files` folder. Here is [an example](project_files/wt_seq.csv). The file should contain **exactly** the two following columns:
- Mutated_seq: all possible values for the Mutated_seq flag from the layout
- WT_seq: corresponding WT DNA sequence, assuming the first three bases constitute the first mutated codon

Alternatively, you can provide the same dataframe with an additional 'nt_seq' column in which you have generated all expected sequences. In this case, the file should be named `expected_mutants.csv` file.

### Experimental design mode

This is a fairly recent feature of gyōza and will be even more simplified in future versions, but is currently **essential**. Choose between 3 possible modes:
- codon: Automatically generate all expected mutants based on the codon mode
- provided: You provide the list of expected mutants or the dataframe of barcode-variant associations (see below)
- random: Random mutagenesis - mutants observed in the sequencing data are directly annotated and filtered based on an acceptable number of amino acid changes

### Codon table

To prevent any typing mistake, the genetic code is imported from a [CoCoPUTs](https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs) table (which also features codon frequencies, although the workflow does not make use of this). [The one provided](project_files/ScerevisiaeTAXID559292_Cocoputs_codon_table.csv) corresponds to *Saccharomyces cerevisiae* TAXID 559292. Please edit the main [config file](config.yaml) if you ever need to specify a different genetic code. Any csv-formatted file with at least two columns ("codon" and "aminoacid") should do.

### Normalization with the number of cellular generations

This normalization is **optional**. Please set the corresponding parameter to True or False in the main config file (see section below). In any case, a csv-formatted file will be **automatically generated** the first time the workflow is run (even if it is a dry run). Again, if normalization is set to True in the config, you will be prompted to edit the file to add the number of cellular generations for each condition in the column 'Nb_gen'. Once the file is edited, re-run the workflow.

### Codon mode

Only needed when the design is set to 'codon'. Please edit the corresponding config entry to specify the type of degenerate codons you introduced at each position in the specified loci. Currently supported are: "NNN" (default value) or "NNK" for single mutants, "NNN x NNN" or "NNK x NNK" for double mutants (including the corresponding single mutants). This is only used to generate the expected sequences for non-barcoded designs. In general, it is preferred to select the 'random' mode and filter out mutants with more than 1 or 2 amino acid changes. This will skip the generation of expected mutants, looking only at the sequencing dataset, resulting in a much faster and less greedy workflow (in terms of computational resources). This functionality will stay for now because in some cases it allows for faster troubleshoot (for example if a library as no coverage at specific positions).

### Barcoded design

To specify a barcoded design:
- Set the "rc_level" parameter of the config to "barcode" (specifies the level at which gyōza should attribute read counts)
- Enter your "barcode attributes", e.g. ['barcode'] or ['bc_index','barcode']
- Provide a datraframe containing all associations of barcodes (column 'barcode') and nucleotide sequences (column 'nt_seq'). This dataframe should respect exactly the same format as `wt_seq.csv` (see section "WT DNA sequences" above), which means in total you should have at least four columns (Mutated_seq, WT_seq, nt_seq, barcode). Additional columns correspond to the other barcode attributes (for example bc_index if we follow the example stated above). The file should be named `expected_mutants.csv` and be located in the `config/project_files` folder.

> [!IMPORTANT]
> 
> For a barcoded design, it is essential to provide the dataframe of barcode-variants association. Since this file already contains all expected sequences, you don't need to provide the `wt_seq.csv` file as well.

Upon completion of the workflow, barcode-level information will be preserved in 'results/df/all_scores.csv', while fitness values will be calculated by aggregating on high-confidence variants (which does not preserve neither barcode-level nor codon-level information).

### Random mutagenesis

To specify a random mutagenesis experimental design:
- Make sure the design entry = 'random'
- Enter an acceptable number of amino acid changes (e.g. 2 = only mutants with up to 2 amino acid changes will be kept)

## Main config file

The main config file is located [here](config.yaml). Please make sure to:
* select the samples to be processed (or leave 'all' if you want to process all samples)
* list your sample attributes
* replace all parameter values with the ones adapted for your project. Note: a first pass might be necessary to establish what would be a good **read count threshold**. Feel free to adjust it and re-run the workflow (if nothing else has changed, only the last steps should run again). This parameter is important because the "avg_scores" dataframe is built only upon "high confidence" variants, i.e. variants with a read count above the set threshold in all T0 replicates.
* set the "perform qc" parameter to True if you want to analyze your raw FASTQ with FastQC (and generate a MultiQC report)
* set the "normalize with gen" parameter to True if you want to normalize with the number of cellular generations
* set the "generate report" parameter to True if you want the HTML report to be automatically generated upon full completion of the workflow
* edit all directory/file paths if necessary

## Note on validation

Currently, all the following files are validated against a YAML schema to help spot formatting issues (misspelled column headers, missing mandatory properties, improper format, etc.): main config file, sample layout, file with WT DNA sequences, codon table, file with the number of cellular generations.

## Profiles for execution

> [!IMPORTANT]
> 
> By default, the simple command line `snakemake` will run gyōza with [the default profile](../profiles/default/config.v8+.yaml) = local execution
> 
> To switch to the SLURM executor, edit [the slurm profile](../profiles/slurm/config.v8+.yaml), including to indicate your email adress

Flags added to the snakemake command line will supersede the values specified in either profile.

> [!WARNING]
> 
> By default, an email will be sent every time a job fails. This is useful to catch TIMEOUT and MEM_OUT errors, but we recommend automatically redirecting emails to prevent inbox overflow.
