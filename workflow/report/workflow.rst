gyōza analyzes deep-mutational scanning (DMS) data to convert allele frequencies of all single mutants of specified loci measured at several time points into functional impact scores.
The main steps are :

1. Trimming with `Cutadapt <http://cutadapt.readthedocs.io>`_ using constant regions
2. Merging with `PANDAseq <https://github.com/neufeld/pandaseq>`_
3. Aggregating with `VSEARCH <https://github.com/torognes/vsearch>`_
4. Comparing non-singleton variants to expected variants
5. Transform read counts into functional impact scores