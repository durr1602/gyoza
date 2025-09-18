[![Conda](https://img.shields.io/badge/conda-≥24.9.1-brightgreen.svg)](https://github.com/conda/conda)
[![Snakemake](https://img.shields.io/badge/snakemake-≥9.4.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/durr1602/gyoza/actions/workflows/main.yml/badge.svg?)](https://github.com/durr1602/gyoza/actions/workflows/main.yml/badge.svg?)
[![DOI](https://zenodo.org/badge/858202238.svg?branch=main&kill_cache=1)](https://zenodo.org/badge/latestdoi/858202238)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

# gyōza: a Snakemake-based workflow to analyze DMS data
<p align="left"><img src="./gyoza.png" width="200"></p>

## Usage

> [!IMPORTANT]
> 
> 1. Read the [documentation](https://durr1602-gyoza.readthedocs.io) and **edit [the main config](config/config.yaml)**
> 2. Choose the type of execution (local or SLURM)

### Check pipeline
3. Perform a dry run using: `snakemake -n`

### Run pipeline
4. Run the workflow either locally: `snakemake` or send to SLURM: `snakemake --profile profiles/slurm`.

## Citation
If you use gyōza, please cite: [doi:10.1093/genetics/iyaf199](https://doi.org/10.1093/genetics/iyaf199).
