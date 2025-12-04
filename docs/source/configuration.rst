Configuration of gyōza
======================

gyōza provides a toy dataset to test installation. Simply run ``snakemake``.

To analyze your data:

1. Edit the main config file: ``config/config.yaml`` with the text editor of your choice
2. Provide the necessary :ref:`project-specific files <project-specific-files>`.

.. _config-file:

Config file
-----------

Below is the config file for the toy dataset.
Edit accordingly to analyze your own data.

.. literalinclude:: ../../config/config.yaml
   :language: yaml
   :caption: Config file for the toy dataset

Make sure you select well your read count threshold (``rc_threshold``),
because this will have an impact on how mutational effects are estimated.

.. _design:

Experimental design
-------------------

In the main config file, under ``project``, choose between 4 possible designs:

- ``codon``: Automatically generate all expected mutants based on the codon mode
- ``provided``: You provide the list of expected mutants
- ``barcode``: You provide the dataframe of barcode-variant associations
- ``random``: Random mutagenesis - mutants observed in the sequencing data are directly
  annotated and filtered based on an acceptable number of amino acid changes (please
  specify this number in the config under ``random``)

Once you’ve selected your design, read carefully what follows to know which files are
needed and which config entries are necessary to edit.

.. _project-specific-files:

Project-specific files
----------------------

Unless stated otherwise, please place your project-specific files in the project folder.
The default path is ``config/project_files/`` but you can edit this in the config under
``project``.

.. _sequencing-data:

Sequencing data
~~~~~~~~~~~~~~~

Please provide the raw reads (forward and, optionally, reverse) of your DMS sequencing
data in the ``config/reads`` folder (or specify a different path in the config under
``reads``). The file names should be featured in the :ref:`layout <layout>`. In the config,
specify if you have provided paired-end reads or not (same type for all samples).

.. tip::

    For paired-end reads, if the mutated DNA sequence can be found in both forward and
    reverse reads, it might be more suitable to analyze only R1 (this is particularly
    recommended for very short mutated sequences). Simply disable the ``paired`` parameter
    under ``reads`` in the config. In the layout, provide linked adapters for
    ``N_forward`` as indicated in the next section.

.. _layout:

Layout
~~~~~~

Please provide a csv-formatted layout of your samples. The file should be named
``layout.csv`` and be located in the project folder. Here is an example:

.. csv-table:: Sample layout
    :header-rows: 1
    :stub-columns: 1
    :file: ../../config/project_files/layout.csv

The file should contain the following columns:

- ``Sample_name``: the unique identifier for each of your samples.
- ``R1``: base name of the fastq file for forward (R1) reads (can be gzipped), including
  extension
- ``R2``: base name of the fastq file for reverse (R2) reads (can be gzipped), including
  extension. Leave empty if you provide single-end sequencing data.
- ``N_forward``: the 5’-3’ DNA sequence corresponding to the fixed region upstream of the
  mutated sequence or anything that can be used as ``-g`` flag with ``cutadapt`` (including
  complex patterns such as ``‘NNATG;optional…ATG’``, in which case do not forget the single
  quotes). For single-end sequencing data (or if you've intentionally specified to analyze
  only R1, as mentioned in the tip above), please specify both constant sequences
  upstream and downstream (on the same strand) separated by ``…``, e.g. ``AAAAGCTG…GCGCTAAAT``
  (no need for single quotes). Note that the use of ``…`` instructs both 5' and 3' trimming
  with linked adapters, which may interfere with the merging step downstream. Therefore this
  option is usually reserved when only R1 is trimmed (in which case, there is no merging).
- ``N_reverse``: the 5’-3’ DNA sequence corresponding to the fixed region 5’ of the mutated
  sequence on the reverse strand or anything that can be used as ``-G`` flag with cutadapt
  (same requirements as above). Leave empty if you provide single-end sequencing data or
  if you're only analyzing R1.
- ``Mutated_seq``: the unique identifier for the mutated DNA sequence, should be the same
  for all samples in which the same sequence was mutated
- ``Pos_start``: starting position in the protein sequence. If you’ve mutated several
  regions/fragments in a coding gene, this position should refer to the full-length
  protein sequence
- ``Replicate``: e.g. ``R1``
- ``Timepoint``: ``T0``, ``T1``, ``T2``, etc. Please provide at least one T0 sample,
  other time points are optional. **T0 samples should not be duplicated for each group**.
- ``Analyze``: ``y`` (or a different truthy value) to process the sample. Leave empty or enter
  non-truthy value to exclude from analysis. Column ignored when ``process_all_samples``
  is enabled in the config.
- ``Report``: ``y`` (or a different truthy value) to include the sample in the HTML report.
  Leave empty or enter non-truthy value to exclude from the report. Samples marked for
  reporting are automatically analyzed. Column ignored when ``process_all_samples``
  is enabled in the config.

Finally, additional columns can be added by the user to specify what makes this sample
unique (other than ``Replicate`` and ``Timepoint``).

Several or all of these additional columns should correspond to what we refer to as:

- sample attributes: what makes an unscreened library unique (minimally ``Mutated_seq``)
- screening attributes: what makes a screened library unique

Please specify both by listing the corresponding columns in the config under ``project``.

Sample attributes may include ``Mutated_seq`` or a combination of attributes that recapitulate
``Mutated_seq`` (as illustrated by the toy dataset). Screening attributes on the other hand
typically include the selective pressure (e.g. ``Drug``). In the sample layout, for the columns
that corresponding to screening attributes, please leave empty for all T0 samples.

In summary, a “sample” is any unique combination of ``Replicate`` + ``Timepoint`` + ``sample
attributes`` + ``screening attributes`` and should be associated to 1 or 2 fastq files,
for the forward and reverse reads, respectively.

.. _codon-table:

Genetic code
~~~~~~~~~~~~

To prevent any typing mistake, the genetic code is imported from a `CoCoPUTs
<https://dnahive.fda.gov/dna.cgi?cmd=codon_usage&id=537&mode=cocoputs>`__ table (which
also features codon frequencies, although the workflow does not make use of this). The
one provided corresponds to *Saccharomyces cerevisiae* TAXID 559292. Feel free to
replace it if you ever need to specify a different genetic code. Any csv-formatted file
with at least two columns (``codon`` and ``aminoacid``) should do. The file should be named
``codon_table.csv`` and be placed in the project folder.

.. _wt-file:

WT DNA sequences
~~~~~~~~~~~~~~~~

If you’ve selected a ``codon`` or ``random`` :ref:`design <design>`, please provide a csv-formatted list of
WT DNA sequences. The file should be named ``wt_seq.csv`` and be located in the project
folder. Here is an example:

.. csv-table:: WT
    :header-rows: 1
    :widths: 1,1,1
    :file: ../../config/project_files/wt_seq.csv

The file should contain the following columns:

- ``Mutated_seq``: all possible values for the Mutated_seq flag from the layout (no
  duplicates!)
- ``WT_seq``: corresponding WT DNA sequence, assuming the first three bases constitute the
  first mutated codon (no duplicates!)

For ``codon`` designs, please add a third column:

- ``codon_mode``: type of degenerate codons you introduced at each position in the locus
  that features on the same row (choose between the currently supported options: ``NNN``,
  ``NNK``, ``NNN x NNN`` or ``NNK x NNK``). This is used to generate the expected sequences.

.. _exp-mut-file:

List of expected mutants
~~~~~~~~~~~~~~~~~~~~~~~~

If you’ve selected the ``provided`` :ref:`design <design>`, please provide 1 compressed dataframe for each
mutated locus, listing all expected sequences. The files should be named
``{Mutated_seq}.csv.gz`` (where ``{Mutated_seq}`` is replaced with the actual label,
e.g. ``Fragment1``) and be located in a subfolder ``expected_mut/``, placed in the project
folder. Each file should contain at least three columns:

- ``Mutated_seq``: a single value per file (out of those listed in the ``Mutated_seq`` column of
  :ref:`the sample layout <layout>`)
- ``WT_seq``: corresponding WT DNA sequence (single value per file), assuming the first
  three bases constitute the first mutated codon
- ``nt_seq``: expected sequences (one per row)

For barcoded :ref:`designs <design>`, please provide the same files with at least one additional column:

- ``barcode``: barcode sequences (one per row, no duplicates!)

Additional columns can be further added to label barcodes with ``barcode attributes`` (for
example, if you want to label each barcode with unique indexes or identifiers). These
barcode attributes can be specified in the config under ``barcode``. Upon completion of
the workflow, barcode-level information will be preserved in
``results/df/all_scores.csv``, while fitness values will be calculated by aggregating on
high-confidence variants (which does not preserve neither barcode-level nor codon-level
information).

Please do not include any index column in the file. Whether the file is actually compressed
or not is not relevant, as long as the file is properly named (ending with ``.csv.gz``).

.. _norm-gen:

Normalization with the number of cellular generations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This normalization is **optional**. Please set the corresponding parameter to ``True`` or
``False`` in the config. In any case, a csv-formatted template named ``nbgen.csv`` will be
**automatically generated** in the project folder the first time the workflow is run
(even if it is a dry run). If normalization is set to ``True`` in the config, you will be
prompted to edit the file to add the number of cellular generations for each condition
(based on current sample selection) in the column ``Nb_gen``. The value entered should
correspond to the number of cellular generations between ``T0`` and the time point on the
matching row. Once the file is edited, re-run the workflow.

.. tip::

    Even if you don’t opt in for this normalization, the generated template can be
    useful to spot any error related to setting up the sample layout.

    Check that the expected groups are listed based on your current selection, with the
    appropriate values for each of your sample and screening attributes.

Final checklist for the main config file
----------------------------------------

Go over your main config file one last time and check the following:

.. |check| unicode:: ☑

- |check| list your ``sample attributes`` and ``screening attributes``
- |check| replace all parameter values with the ones adapted for your project. Note: a
  first pass might be necessary to establish what would be a good **read count
  threshold** (specified under ``reads``). Feel free to adjust it and re-run the workflow
  (if nothing else has changed, only the last steps should run again). This parameter is
  important because the ``avg_scores`` dataframe is built only upon “high confidence”
  variants, i.e. variants with a read count above the set threshold in all T0
  replicates.
- |check| set the ``process_all_samples`` to ``True`` if you want to include all samples for
  processing and reporting
- |check| set the ``perform_qc`` parameter to ``True`` if you want to analyze your raw FASTQ
  with Fastp (and generate a MultiQC report)
- |check| set the ``process_frequencies`` to ``True`` if you want to convert allele frequencies
  into functional impact scores (``False`` if you simply want heatmaps of raw read counts and
  frequencies, e.g. to assess diversity in unscreened libraries)
- |check| set the ``normalize_with_gen`` parameter to ``True`` if you want to normalize with
  the number of cellular generations (only valid if you opted in for processing frequencies)
- |check| edit the directory paths to :ref:`project-specific files <project-specific-files>`
  and :ref:`reads <sequencing-data>` if necessary.

Note on validation
------------------

Currently, all the following files are validated against a YAML-formatted JSON schema to
help spot formatting issues (misspelled column headers, missing mandatory properties,
improper format, etc.):

- :ref:`main config file <config-file>`
- :ref:`sample layout <layout>`
- :ref:`file with WT DNA sequences <wt-file>`
- :ref:`files with expected sequences <exp-mut-file>`
- :ref:`codon table <codon-table>`
- :ref:`file with the number of cellular generations <norm-gen>`

These validations **do not cover** the following:

- if CSV files are actually formatted as TSV
- if there's a problem with specific row(s) of a CSV file

In some cases, the error message displayed might still help you troubleshoot.

Profiles for execution
----------------------

.. important::

    By default, the simple command line ``snakemake`` will run gyōza with the default
    profile (``profiles/default/config.v8+.yaml``) = local execution

    To switch to the SLURM executor, edit the slurm profile
    (``profiles/slurm/config.v8+.yaml``), including to indicate your email address

Flags added to the ``snakemake`` command line will supersede the values specified in either
profile.

.. warning::

    By default, an email will be sent every time a job fails. This is useful to catch
    ``TIMEOUT`` and ``MEM_OUT`` errors, but we recommend automatically redirecting emails to
    prevent inbox overflow.
