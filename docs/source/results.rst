Exploring the results
=====================

Main report
-----------

gyōza automatically generates an HTML report, regardless of the state of
completion of the workflow, which means it is built from the files
generated over the last run (and/or previous runs if it matches the config).

Depending on the number of files included, the report may be compressed
but should be located at the root of the ``results/`` folder.

Depending on the config and what was successfully generated, the following
sections may be included in the report.

Main page
~~~~~~~~~

On the main page of the report, one should find:

- Metadata for the workflow (information related to the config)
- The rule graph describing the workflow

Following hyperlinks in the metadata will lead to the documentation of
external tools used by gyōza.

Clicking on the nodes in the rule graph will display the **rule definitions**
in the side panel on the left. For each rule, this includes:

- Input
- Output
- Software (packages and versions)
- Source code

Side panel
~~~~~~~~~~

.. |home-icon| raw:: html

   <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" style="width: 1.5em; height: 1.5em;">
     <path stroke-linecap="round" stroke-linejoin="round" d="m2.25 12 8.954-8.955c.44-.439 1.152-.439 1.591 0L21.75 12M4.5 9.75v10.125c0 .621.504 1.125 1.125 1.125H9.75v-4.875c0-.621.504-1.125 1.125-1.125h2.25c.621 0 1.125.504 1.125 1.125V21h4.125c.621 0 1.125-.504 1.125-1.125V9.75M8.25 21h8.25" />
   </svg>

.. |eye-icon| raw:: html

   <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" style="width: 1.5em; height: 1.5em; vertical-align: middle;">
     <path stroke-linecap="round" stroke-linejoin="round" d="M2.036 12.322a1.012 1.012 0 0 1 0-.639C3.423 7.51 7.36 4.5 12 4.5c4.638 0 8.573 3.007 9.963 7.178.07.207.07.431 0 .639C20.577 16.49 16.64 19.5 12 19.5c-4.638 0-8.573-3.007-9.963-7.178Z" />
     <path stroke-linecap="round" stroke-linejoin="round" d="M15 12a3 3 0 1 1-6 0 3 3 0 0 1 6 0Z" />
   </svg>

.. |info-icon| raw:: html

   <svg xmlns="http://www.w3.org/2000/svg" fill="none" viewBox="0 0 24 24" stroke-width="1.5" stroke="currentColor" style="width: 1.5em; height: 1.5em; vertical-align: middle;">
   <path stroke-linecap="round" stroke-linejoin="round" d="m11.25 11.25.041-.02a.75.75 0 0 1 1.063.852l-.708 2.836a.75.75 0 0 0 1.063.853l.041-.021M21 12a9 9 0 1 1-18 0 9 9 0 0 1 18 0Zm-9-3.75h.008v.008H12V8.25Z" />
   </svg>

In the side panel, regardless of what opens in the main view, the
|home-icon| button brings you back home.

The "home" page for the side panel corresponds to "Workflow" under
**General**. Statistics include run time and date of jobs for each rule.
A search bar is available at the very top from the "home" page.
Results are accessible under **Results** and can include any of the
following sections.

.. important::

    Whenever you see |eye-icon| |info-icon|, clicking on |eye-icon| will
    refresh the main view by loading the corresponding plot, while clicking
    on |info-icon| will display more information, typically the caption for
    the plot.

Quality control
~~~~~~~~~~~~~~~

The "Quality control" section currently corresponds to a single entry,
which is "Interactive QC report". Clicking on |eye-icon| should open the
MultiQC report. **If it doesn't**, you can still open the MultiQC report
separately (``results/0_qc/multiqc.html``).

Read filtering
~~~~~~~~~~~~~~

The "Read filtering" section contains up to two entries:

- Aggregated: contains plots with several samples on each plots

  - Summary of filtered reads
  - Read counts of unexpected variants  

- Heatmaps of raw read counts

Read processing
~~~~~~~~~~~~~~~

The "Read processing" section appears only when you've enabled
``process_read_counts`` in the config. It contains:

- one entry per unique combination of :ref:`sample attributes <layout>`:

  - Raw read count per variant
  - Overlap across time points and replicates

- Aggregated: contains plots with pooled groups of samples

  - Distribution of allele frequencies

Functional impact
~~~~~~~~~~~~~~~~~

The "Functional impact" sections appears only when you've enabled
``process_read_counts`` in the config. It contains:

- Correlation between time points
- Aggregated
  
  - Distribution of functional impact scores
  - Correlation between replicates (1/2): heatmap of Spearman correlation
    coefficients
  - Correlation between replicates (2/2): scatterplot of pairwise
    comparisons
  - Functional impact over time

- Heatmaps of functional impact

Generated files
---------------

Some files are two heavy to be kept and are therefore automatically removed.
You can disable this option by running the workflow with the ``--no-temp``
flag.

All other non-log files are located in the ``results/`` folder.

Fastp and MultiQC reports can be found in ``results/0_qc``.

Dataframes can be found in the ``results/df``, including:

- ``all_scores.csv`` which contains all nonaggregated scores
- ``avg_scores.csv`` which contains fitness and error values for
  high-confidence variants

To obtain this last dataframe, gyōza calculates the median functional
impact score for each amino acid sequence (from high-confidence variants
only), then calculates the median and error across replicates. The lower
and upper error values (``lower_err`` and ``upper_err``, respectively)
are obtained by subtracting the 2.5th or 97.5th percentile.

.. note::

    For ``barcode`` designs, it is often useful to show what happens
    when we sample ``n`` barcodes per mutation. In any case, these
    additional downstream steps should be manually performed from
    ``all_scores.csv``, which preserves barcode-level information.

Graphs can be found in the ``results/graphs`` folder.