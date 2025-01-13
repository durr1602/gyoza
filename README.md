# Snakemake workflow: Analysis of DMS data

[![Conda](https://img.shields.io/badge/conda-≥24.9.1-brightgreen.svg)](https://github.com/conda/conda)
[![Snakemake](https://img.shields.io/badge/snakemake-≥8.23.2-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/durr1602/DMS_analysis_snakemake/workflows/Tests/badge.svg?branch=main)](https://github.com/durr1602/DMS_analysis_snakemake/actions?query=branch%3Amain+workflow%3ATests)


A Snakemake workflow to analyze demultiplexed sequencing data of a DMS experiment, produce diagnostics plots and generate a dataframe of selection coefficients.

## Installation

0. If you use Windows, the following steps need to be run in WSL2 ([WSL installation instructions](https://learn.microsoft.com/en-us/windows/wsl/install)). If you already have WSL2 or if you use Linux or MacOS X, the following steps need `conda`. If you don't already have it, we recommend installing Miniforge by following the instructions listed in the "Step 1" section of [this tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-1-installing-miniforge). Continue by running the command lines listed here (clone repo, ...) in the same terminal.

1. Clone this repository:
```
git clone https://github.com/durr1602/DMS_analysis_snakemake.git
```
2. Move to the corresponding directory:
```
cd DMS_analysis_snakemake
```
3. If you correctly followed instructions at step 0, you should be able to use `conda`, first update to a recent version (>=24.7.1, ideally even more recent such as >=24.9.1). You might run into some issues, so one workaround is to update using `mamba`:
```
mamba update conda
```

4. Install snakemake in a virtual environment (the "solving environment" step can take some time, normally not longer than a few minutes). Again, you could run into some issues if you have multiple instances of `conda` on your system (e.g. anaconda and miniforge), in which case you can uninstall (for example) anaconda by running `rm -rf anaconda3`.
```
conda env create --name=DMS-snake --file=env.yml
conda activate DMS-snake
```

## Usage

### Prepare files and edit config
1. **IMPORTANT**: Read the [config documentation](config/README.md) and **edit the main config**. If you plan on sending the pipeline to SLURM, make sure you also **edit the technical config file**.

### Check pipeline
2. (recommended) Perform a dry run using: `snakemake --use-conda -n`

This step is strongly recommended. It will make sure the prepared workflow does not contain any error and will display the rules (steps) that need to be run in order to reach the specified target(s) (default value for the target is the dataframe of selection coefficients, which is produced during the very last step of the workflow). For both the dry run and the actual run, you can decide to run the workflow only until a certain file is generated or rule is completed, using the `--until` flag in the snakemake command line, for example: `snakemake -n --until stats`

### Run pipeline
3. Running the workflow with `conda`

    a) Locally: `snakemake --cores 4 --use-conda` (recommended only for small steps or to run the workflow on the provided example dataset, with the `--cores` flag indicating the max number of CPUs to use in parallel - can be adapted depending on the resources available on your machine).
    
    b) **or** send to SLURM (1 job per rule per sample): `snakemake --profile profile --sdm conda` (all parameters are specified in the [tech config file](profile/config.v8+.yaml), jobs wait in the queue until the resources are allocated. For example, if you're allowed 40 CPUs, only 4 jobs at 10 CPUs each will be able to run at once. Once those jobs are completed, the next ones in the queue will automatically start.

4. Running the workflow inside a container (needs `singularity`/`apptainer`). In this case, the container will first be created, then conda envs will be created for each rule inside the container. This option is meant to be used on a system where you want to isolate the (many) files installed by `conda`. This option is **not** suited for local execution. The command to send to SLURM (again, 1 job per rule per sample) is: `snakemake --profile profile --sdm conda apptainer` (comments about resource allocation in step 3b apply).

Fore more info on cluster execution: read the doc on [smk-cluster-generic plugin](https://github.com/jdblischak/smk-simple-slurm/tree/main)

**Important** If snakemake is launched directly from the command line, the process will be output to the terminal. Exiting with `<Ctrl+C>` is currently interpreted (as specified in the [tech config file](profile/config.v8+.yaml)) as cancelling all submitted jobs (`scancel`). Exiting during a local execution will **also** abort the workflow. For a small run, this should not be too inconvenient, and you can still open a new terminal to get a prompt. For workflows with a longer runtime, one might want to use `tmux` as described below.

To launch the pipeline and ensure that it continues to run in the background even when the terminal is closed, one should use [`tmux`](https://github.com/tmux/tmux/wiki/Getting-Started). Please make sure the tool is installed first (already installed on some servers). Then, follow the steps:
1. Type `tmux new -s snakes` to launch a new tmux session
2. Activate the conda env with `mamba activate DMS-snake` or `conda activate DMS-snake`
3. Navigate to the Snakefile directory and launch the pipeline with `snakemake --profile profile`
4. To close (detach) the session, type `<Ctrl+b>`, then `<d>`. You should see the message: `[detached (from session snakes)]`
5. To reconnect (attach) to the session, for example from a different machine: `tmux attach -t snakes`. You can also see existing sessions with `tmux ls`.
6. To close the session when everything is finished, type `<Ctrl+b>`, then `<:>`, then `kill-session` and finally `<Enter>`.

### Generate HTML report

Once the workflow has run, you can generate a comprehensive HTML report, complete with underlying code and embedded plots:
```
snakemake --report report.html
```

### Edit pipeline
One can manually edit the [Snakefile](workflow/Snakefile) and/or the rules (.smk files in rules folder) to edit the main steps of the pipeline. This should not be required to run the standard pipeline and should be done only when the workflow itself needs to be modified.
    
**Editing template jupyter notebooks** is tricky to do manually because the paths and kernel are not shared between platforms. Thankfully, there is a snakemake command that allows interactive editing of any template notebook, using any output file (from the notebook) as argument. The following example will generate URLs to open `jupyter`, in which we can edit the process_read_counts notebook that outputs the upset_plot.svg file, as specified in the corresponding .smk file.

```
snakemake --use-conda --cores 1 --edit-notebook ../results/graphs/upset_plot.svg
```

**Careful**, if you are running `snakemake` on a server, you might need to open a SSH tunnel between your local machine and the server by running the following command from a local terminal (should not be necessary when running locally on your machine):
```  
ssh -L 8888:localhost:8888 <USER>@<ADRESS>
```
(adapt port if necessary, 8888 or 8889, should match what is featured in the URL generated by snakemake/jupyter)

The command simply opens a terminal on the server, but now you can copy-paste the URL in your local browser.
    
You can then open the notebook, run it (kernel and paths taken care of) and save it once the modifications have been done. Then click on "Close notebook" and finally "Shut down" on the main jupyter dashboard. The quitting and saving should be displayed on the initial server's terminal (the one from which the `snakemake` command was run). You'll also notice that the size of the template notebook file could be smaller, because outputs are automatically flushed out. To retrieve a notebook with outputs (for future HTML export for example), locate the notebook in the appropriate folder once the pipeline has run (path specified in the log directive of the .smk file).

## Issues / work in progress

Known (potential) issues:
* Tests performed by github actions are currently not working with the latest `snakemake` version (no `conda`) but should be soon

Work in progress:
* Portability / installation on DRAC servers
