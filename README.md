[![Conda](https://img.shields.io/badge/conda-≥24.9.1-brightgreen.svg)](https://github.com/conda/conda)
[![Snakemake](https://img.shields.io/badge/snakemake-≥9.4.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/durr1602/gyoza/actions/workflows/main.yml/badge.svg?)](https://github.com/durr1602/gyoza/actions/workflows/main.yml/badge.svg?)
[![DOI](https://zenodo.org/badge/858202238.svg?branch=main&kill_cache=1)](https://zenodo.org/badge/latestdoi/858202238)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# gyōza: a Snakemake-based workflow to analyze DMS data
<p align="left"><img src="./gyoza.png" width="200"></p>

## Installation

> [!TIP]
> 
> Follow [the full installation instructions](fulldoc/README.md).

## Usage

### Prepare files and edit config

> [!IMPORTANT]
> 
> 1. Read the [config documentation](config/README.md) and **edit [the main config](config/config.yaml)**.
> 2. Choose the type of execution (local or SLURM), edit [the slurm profile](profiles/slurm/config.v8+.yaml) if needed.

### Check pipeline
2. Perform a dry run using: `snakemake -n`

This step is strongly recommended. It will make sure the prepared workflow does not contain any error and will display the rules (steps) that need to be run in order to run the workflow (built dynamically based on the config). If you're running the workflow for the first time and you toggled in normalization with growth data, you should see a warning prompting you to edit the generated template file (for more details, go back to step 1).

### Run pipeline
3. Run the workflow either locally: `snakemake` or send to SLURM: `snakemake --profile profiles/slurm`.
