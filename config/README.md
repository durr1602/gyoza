# Configuration of gyōza

gyōza provides a toy dataset to test installation. Simply run `snakemake`.

To analyze your data, edit the main [config file](config.yaml). and provide the necessary project-specific files.

## Experimental design

In the main config file, choose between 4 possible designs:
- codon: Automatically generate all expected mutants based on the codon mode
- provided: You provide the list of expected mutants
- barcode: You provide the dataframe of barcode-variant associations
- random: Random mutagenesis - mutants observed in the sequencing data are directly annotated and filtered based on an acceptable number of amino acid changes

Once you've selected your design, read carefully what follows to know which files are needed and which config entries are necessary to edit.

## Project-specific files

### Sequencing data

Please provide the raw reads (forward and, optionally, reverse) of your DMS sequencing data in the `config/reads` folder (or specify a different path in the main config file). The file names should be featured in the layout (see section below). Don't forget to specify in the config if you have provided paired-end reads or not.

### Layout

Please provide a csv-formatted layout of your samples. The file should be named `layout.csv` and be located in the `config/project_files` folder. Here is [an example](project_files/layout.csv). The file should contain the following columns:
- Sample_name: the unique identifier for each of your samples.
- R1: base name of the fastq file for forward (R1) reads (can be gzipped), including extension
- R2: base name of the fastq file for reverse (R2) reads (can be gzipped), including extension. Leave empty if you provide single-end sequencing data.
- N_forward: the 5'-3' DNA sequence corresponding to the fixed region upstream of the mutated sequence or anything that can be used as -g flag with cutadapt (including complex patterns such as 'NNATG;optional...ATG', in which case do not forget the single quotes). For single-end sequencing data, please specify both constant sequences upstream and downstream (on the same strand) separated by '...', e.g. AAAAGCTG...GCGCTAAAT (no need for single quotes)
- N_reverse: the 5'-3' DNA sequence corresponding to the fixed region 5' of the mutated sequence on the reverse strand or anything that can be used as -G flag with cutadapt (same requirements as above). Leave empty if you provide single-end sequencing data.
- Mutated_seq: the unique identifier for the mutated DNA sequence, should be the same for all samples in which the same sequence was mutated
- Pos_start: starting position in the protein sequence. If you've mutated several regions/fragments in a coding gene, this position should refer to the full-length protein sequence
- Replicate: e.g. "R1"
- Timepoint: "T0", "T1", "T2", etc. Please provide at least 1 "T0" sample per group, other time points are optional.
- Analyze: "y" (or a different truthy value) to process the sample. Leave empty or enter non-truthy value to exclude from analysis. Corresponding T0 samples and matching replicates are automatically rescued, regardless of the selection. In other words, you can select a single replicate for each group you want to analyze.
- Report: "y" (or a different truthy value) to include the sample in the HTML report. Leave empty or enter non-truthy value to exclude from the report. Samples marked for reporting are rescued as describe above and will be automatically analyzed.

Finally, additional columns can be added by the user to specify what makes this sample unique. These are referred to as "sample attributes" and could correspond to the genetic background, the fragment/region of the gene if it applies (in which case, sample attributes can overlap with Mutated_seq, as is the case for the toy dataset), the drug used for selection, etc. In summary, a "sample" is any unique combination of Mutated_seq + Replicate + Timepoint + sample attributes and should be associated to 1 or 2 fastq files, for the forward and reverse reads, respectively.

### Codon table

To prevent any typing mistake, the genetic code is imported from a [CoCoPUTs](https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs) table (which also features codon frequencies, although the workflow does not make use of this). [The one provided](project_files/ScerevisiaeTAXID559292_Cocoputs_codon_table.csv) corresponds to *Saccharomyces cerevisiae* TAXID 559292. Please edit the main [config file](config.yaml) if you ever need to specify a different genetic code. Any csv-formatted file with at least two columns ("codon" and "aminoacid") should do.

### WT DNA sequences

If you've selected a 'codon' or 'random' design, please provide a csv-formatted list of WT DNA sequences. The file should be named `wt_seq.csv` and be located in the `config/project_files` folder. Here is [an example](project_files/wt_seq.csv). The file should contain **exactly** the two following columns:
- Mutated_seq: all possible values for the Mutated_seq flag from the layout (no duplicates!)
- WT_seq: corresponding WT DNA sequence, assuming the first three bases constitute the first mutated codon (no duplicates!)

### List of expected mutants

If you've selected the 'provided' design, please provide 1 compressed dataframe for each mutated locus, listing all expected sequences. The files should be named `{Mutated_seq}.csv.gz` (where `{Mutated_seq}` is replaced with the actual label, e.g. Fragment1) and be located in the `config/project_files/expected_mut` folder. Each file should contain at least three columns:
- Mutated_seq: a single value per file (out of those listed in the Mutated_seq column of the sample layout)
- WT_seq: corresponding WT DNA sequence, assuming the first three bases constitute the first mutated codon (no duplicates!)
- nt_seq: expected sequences (one per row)

For barcoded designs, please provide the same files with at least one additional column:
- barcode: barcode sequences (one per row, no duplicates!)

Additional columns can be further added to label barcodes with "barcode attributes" (for example, if you want to label each barcode with unique indexes or identifiers). Upon completion of the workflow, barcode-level information will be preserved in `results/df/all_scores.csv`, while fitness values will be calculated by aggregating on high-confidence variants (which does not preserve neither barcode-level nor codon-level information).

### Normalization with the number of cellular generations

This normalization is **optional**. Please set the corresponding parameter to True or False in the main config file. In any case, a csv-formatted template will be **automatically generated** the first time the workflow is run (even if it is a dry run). Again, if normalization is set to True in the config, you will be prompted to edit the file to add the number of cellular generations for each condition (based on current sample selection) in the column 'Nb_gen'. Once the file is edited, re-run the workflow.

## Conditional config entries

Some entries in the main config may depend on other config entries, e.g. the design you've specified.

### Codon mode

Only needed for 'codon' designs. Please edit the corresponding config entry to specify the type of degenerate codons you introduced at each position in the specified loci. Currently supported are: "NNN" (default value) or "NNK" for single mutants, "NNN x NNN" or "NNK x NNK" for double mutants (including the corresponding single mutants). This is used to generate the expected sequences for non-barcoded designs.

### Barcoded design

Only needed for 'barcode' designs.
- Set the "rc_level" parameter of the config to "barcode" (specifies the level at which gyōza should attribute read counts)
- Enter your "barcode attributes", e.g. ['barcode'] or ['bc_index','barcode']

### Random mutagenesis

Only needed for 'random' designs.
- Enter an acceptable number of amino acid changes (e.g. 2 = only mutants with up to 2 amino acid changes will be kept)

## Final checklist for the main config file

Go over your [main config file](config.yaml) one last time and check the following:
- [ ] list your sample attributes
- [ ] replace all parameter values with the ones adapted for your project. Note: a first pass might be necessary to establish what would be a good **read count threshold**. Feel free to adjust it and re-run the workflow (if nothing else has changed, only the last steps should run again). This parameter is important because the "avg_scores" dataframe is built only upon "high confidence" variants, i.e. variants with a read count above the set threshold in all T0 replicates.
- [ ] set the "perform qc" parameter to True if you want to analyze your raw FASTQ with FastQC (and generate a MultiQC report)
- [ ] set the "process_read_counts" to True if you want to convert read counts to functional impact scores (False if you simply want read counts, e.g. to assess diversity in T0 libraries)
- [ ] set the "normalize with gen" parameter to True if you want to normalize with the number of cellular generations (only valid if you opted in for processing read counts)
- [ ] edit all directory/file paths if necessary

## Note on validation

Currently, all the following files are validated against a YAML schema to help spot formatting issues (misspelled column headers, missing mandatory properties, improper format, etc.): main config file, sample layout, file with WT DNA sequences, files with expected sequences, codon table, file with the number of cellular generations.

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
