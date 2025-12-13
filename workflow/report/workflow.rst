**Welcome to the gyōza report!**

Config
------

``>`` gyōza version: ``{{ snakemake.config["pipeline_version"] }}``

``>`` gyōza analyzed your DMS data generated with a ``{{ snakemake.config["project"]["design"] }}`` design.

``>`` You specified the following sample attributes : ``{{ snakemake.config["project"]["sample_attributes"] }}``

``>`` You specified the following screening attributes : ``{{ snakemake.config["project"]["screening_attributes"] }}``

Read processing
---------------

{% if snakemake.config["reads"]["paired"] %}
``>`` Trimming was performed on paired-end reads with `Cutadapt <http://cutadapt.readthedocs.io>`_ using constant regions.

``>`` Overlapping reads were merged with `BBMerge <https://github.com/bbushnell/BBTools>`_.
{% else %}
``>`` Trimming was performed on single-end reads with `Cutadapt <http://cutadapt.readthedocs.io>`_ using constant regions.
{% endif %}

``>`` Read counts were aggregated with `VSEARCH <https://github.com/torognes/vsearch>`_.

{% if snakemake.config["project"]["design"] == "random" %}
Read filtering
--------------

``>`` Variants with more than {{ snakemake.config["random"]["Nham_aa_max"] }} amino acid changes were labeled as "Unexpected" and discarded.
{% endif %}

``>`` Read counts were converted into allele frequencies.

{% if snakemake.config["process_frequencies"] %}
Estimate mutational effects
---------------------------

``>`` Read counts were converted into functional impact scores{% if snakemake.config["normalize_with_gen"] %}, with growth normalization{% endif %}.

``>`` You specified a read count threshold of {{ snakemake.config["reads"]["rc_threshold"] }}.
{% endif %}

**Thanks for using gyōza!**