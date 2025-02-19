[![Conda](https://img.shields.io/badge/conda-≥24.9.1-brightgreen.svg)](https://github.com/conda/conda)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.23.2-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/durr1602/gyoza/workflows/Tests/badge.svg?branch=main)](https://github.com/durr1602/gyoza/actions?query=branch%3Amain+workflow%3ATests)

# gyōza: a Snakemake-based workflow to analyze DMS data
<p align="left">
  <img src="./gyoza.png" width="200">
</p>

## Installation

If you are not familiar with cloning a repo / conda environments, etc, or if you like reading, please refer to [the full installation instructions](fulldoc/README.md). Otherwise, after cloning the repo, make sure your version of conda is up to date (>=24.7.1 required) and create the base env from the provided requirements file:
```
conda env create --name=gyoza --file=env.yml
conda activate gyoza
```

## Usage

### Prepare files and edit config
1. **IMPORTANT**: Read the [config documentation](config/README.md) and **edit the main config**. If you plan on sending the pipeline to SLURM, make sure you also **edit the technical config file**.

### Check pipeline
2. Perform a dry run using: `snakemake -n`

This step is strongly recommended. It will make sure the prepared workflow does not contain any error and will display the rules (steps) that need to be run in order to reach the specified target(s) (default targets include the dataframe of functional impact scores, which is produced during the very last step of the workflow). If you're running the workflow for the first time and you toggled in normalization with growth data, you should see a warning prompting you to edit the generated template file (for more details, go back to step 1).

### Run pipeline
3. Run the workflow either locally: `snakemake --use-conda` **or** send to SLURM: `snakemake --profile profile`. In the latter case, make sure to first edit the [tech config file](profile/config.v8+.yaml).
