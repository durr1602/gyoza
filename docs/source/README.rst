.. image:: https://img.shields.io/badge/conda-≥24.9.1-brightgreen.svg
    :target: https://github.com/conda/conda
    :alt:

.. image:: https://img.shields.io/badge/snakemake-≥9.4.0-brightgreen.svg
    :target: https://snakemake.github.io
    :alt:

.. image:: https://zenodo.org/badge/858202238.svg?branch=main&kill_cache=1
    :target: https://zenodo.org/badge/latestdoi/858202238
    :alt:

gyōza: a Snakemake-based workflow to analyze DMS data
=====================================================

.. image:: ../../gyoza.png
    :height: 100px
    :alt: gyoza logo
    :align: left

gyōza is a command-line interface to help you analyze time-series deep-mutational
scanning (DMS) data .

After `installation <installation.html>`__, prepare the `configuration
<configuration.html>`__ and `run the workflow <usage.html>`__.

DMS sequencing data are processed with the following existing software:

- `FastQC <https://github.com/s-andrews/FastQC>`__
- `Cutadapt <http://cutadapt.readthedocs.io>`__
- `PANDAseq <https://github.com/neufeld/pandaseq>`__
- `VSEARCH <https://github.com/torognes/vsearch>`__

gyōza was developed by `Romain Durand <mailto:duran2101@gmail.com>`__
|ORCID_icon|
during his postdoc in the `Landry lab <https://landrylab.ibis.ulaval.ca/>`__
at `Université Laval <https://www.ulaval.ca/>`__.

Citation
~~~~~~~~
If you use gyōza, please cite: `doi:10.1093/genetics/iyaf199 <https://doi.org/10.1093/genetics/iyaf199>`__.

.. |ORCID_icon| image:: _static/iD_icon.png
   :target: https://orcid.org/0000-0002-7681-4727
   :alt: ORCID profile
