|Conda| |Snakemake| |DOI|

gyōza: a Snakemake-based workflow to analyze DMS data
=====================================================

.. raw:: html

   <p align="left">

.. raw:: html

   </p>

Installation
------------

   [!TIP]

   Follow `the full installation instructions <fulldoc/README.md>`__.

Usage
-----

Prepare files and edit config
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   [!IMPORTANT]

   1. Read the `config documentation <config/README.md>`__ and **edit**
      `the main config <config/config.yaml>`__.
   2. Choose the type of execution (local or SLURM), edit `the slurm
      profile <profiles/slurm/config.v8+.yaml>`__ if needed.

Check pipeline
~~~~~~~~~~~~~~

2. Perform a dry run using: ``snakemake -n``

This step is strongly recommended. It will make sure the prepared
workflow does not contain any error and will display the rules (steps)
that need to be run in order to run the workflow (built dynamically
based on the config). If you’re running the workflow for the first time
and you toggled in normalization with growth data, you should see a
warning prompting you to edit the generated template file (for more
details, go back to step 1).

Run pipeline
~~~~~~~~~~~~

3. Run the workflow either locally: ``snakemake`` or send to SLURM:
   ``snakemake --profile profiles/slurm``.

.. |Conda| image:: https://img.shields.io/badge/conda-≥24.9.1-brightgreen.svg
   :target: https://github.com/conda/conda
.. |Snakemake| image:: https://img.shields.io/badge/snakemake-≥9.4.0-brightgreen.svg
   :target: https://snakemake.github.io
.. |DOI| image:: https://zenodo.org/badge/858202238.svg?branch=main&kill_cache=1
   :target: https://zenodo.org/badge/latestdoi/858202238
