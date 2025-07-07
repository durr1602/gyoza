# Documentation for gyōza

## Full installation instructions

If you use Windows, the following steps need to be run in WSL2 ([WSL installation instructions](https://learn.microsoft.com/en-us/windows/wsl/install)).

If you already have WSL2 or if you use Linux or MacOS X, the following steps need `conda` (see [Conda documentation](https://conda.io/docs/index.html)).

If you don't already have it, we recommend installing Miniforge by following the instructions listed in the "Step 1" section of [this tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-1-installing-miniforge).

If you already have it, make sure to update it to a recent version (>=24.7.1, ideally even more recent such as >=24.9.1). Here are two helpful command lines to do this: `conda update -n base -c defaults conda --repodata-fn=repodata.json`, or `mamba update conda`. If you already have it or another instance such as Anaconda, you can uninstall (for example) anaconda by running `rm -rf anaconda3`.

1. Create a Python>=3.9 virtual environment to install dependencies. You may need to adapt `python3.11` depending on your version (try `python --version` or `python3 --version` for example, to check the Python version).
```
python3.11 -m venv gyoza_env
source gyoza_env/bin/activate
pip install snakedeploy>=0.11.0 snakemake>=9.4.0 snakemake-wrapper-utils>=0.7.2 pygments>=2.19.1 snakemake-executor-plugin-cluster-generic
```

### Recommended installation of gyōza

2. Deploy gyōza.

The following command line uses [Snakedeploy](https://snakedeploy.readthedocs.io/en/latest/index.html) to create a minimal file tree that will allow you to use gyōza.

It is specific to both the version of gyōza deployed and the DMS project your want to analyze, and should make it easier to create a repository for improved reproducibility. The --tag argument accepts any branch or release version tag
```
snakedeploy deploy-workflow https://github.com/durr1602/gyoza my_gyoza_project --tag main
cd my_gyoza_project
```

### Alternative option: clone the repository

This option provides the most flexibility, since the entire repository is cloned, meaning one can potentially modify the code of gyōza.

2. Clone this repository:
```
git clone https://github.com/durr1602/gyoza.git
cd gyoza
```

## Usage

### Prepare files and edit config

> [!IMPORTANT]
> 
> 1. Read the [config documentation](../config/README.md)
> 2. **Edit [the main config](../config/config.yaml)**
>
> If you choose to execute gyōza locally, you may optionally edit [the default profile](../profiles/slurm/config.v8+.yaml)
> 
> If you choose to execute gyōza with SLURM, don't forget to edit [the slurm profile](../profiles/slurm/config.v8+.yaml)

### (optional) Prepare environments

> [!TIP]
>
> 3. (optional) Create all conda environments using: `snakemake --conda-create-envs-only`.
>
> This step is only required before first use and is always included when running the workflow. Unfortunately, at the time of writing, validations will be run even for these command lines. This means that you need to fully prepare the workflow before creating all envs.

### Check pipeline

> [!IMPORTANT]
> 
> 4. Perform a dry run using: `snakemake -n`
>
> This step is strongly recommended. It will make sure the prepared workflow does not contain any error and will display the rules (steps) that need to be run in order to run the workflow (built dynamically based on the config). If you're running the workflow for the first time and you toggled in normalization with growth data, you should see a warning prompting you to edit the generated template file (for more details, go back to step 1).

### Run pipeline

> [!IMPORTANT]
> 
> 5. Run the workflow either **locally (a)** or **using SLURM (b)**
> a) Local execution: `snakemake` or `snakemake --profile profiles/default`
> b) SLURM (1 job per rule per sample): `snakemake --profile profiles/slurm`. Jobs wait in the queue until the resources are allocated. For example, if you're allowed 40 CPUs, only 4 jobs at 10 CPUs each will be able to run at once. Once those jobs are completed, the next ones in the queue will automatically start. Fore more info on cluster execution: read the doc on [smk-cluster-generic plugin](https://github.com/jdblischak/smk-simple-slurm/tree/main)

### Abort pipeline / exit terminal

If snakemake is launched directly from the command line, the process will be output to the terminal. Exiting with `<Ctrl+C>` is currently interpreted (as specified in the [the slurm profile](../profiles/slurm/config.v8+.yaml)) as cancelling all submitted jobs (`scancel`). Exiting during a local execution will **also** abort the workflow. This means that while the workflow is running, the user cannot get the prompt back.

There are 3 possible options to get the prompt back and/or exit the terminal without aborting the workflow:
1. Open a new tab on your terminal (may require to log into the session again)
2. Use `nohup` (e.g. `nohup snakemake`). Closing the tab will not abort the workflow.
3. Use the terminal multiplexer [`tmux`](https://github.com/tmux/tmux/wiki/Getting-Started). This way you can get the prompt back right away, exit without aborting and even reconnect to the session from a different machine.

Please make sure `tmux` is installed (already installed on some servers). Then, follow the steps:
1. Type `tmux new -s snakes` to launch a new tmux session
2. Activate the conda env with `mamba activate gyoza` or `conda activate gyoza`
3. Navigate to the Snakefile directory and launch the pipeline with `snakemake --profile profiles/slurm`
4. To close (detach) the session, type `<Ctrl+b>`, then `<d>`. You should see the message: `[detached (from session snakes)]`
5. To reconnect (attach) to the session, for example from a different machine: `tmux attach -t snakes`. You can also see existing sessions with `tmux ls`.
6. To close the session when everything is finished, type `<Ctrl+b>`, then `<:>`, then `kill-session` and finally `<Enter>`.

### Useful Snakemake command lines

> [!TIP]
> 
> For both the dry run and the actual run, you can decide to run the workflow only until a certain file is generated or rule is completed, using the `--until` flag in the snakemake command line
>
> For example: `snakemake -n --until stats`

>[!TIP]
>
> Similarly, you can omit rules from the workflow.
>
> For example: `snakemake --omit rule_to_omit`

>[!NOTE]
>
> To account for intermediate runs, the report is built dynamically based on existing files (out of all files expected to be included in the report). If any error were to arise, one can still generate an intermediate report with the dedicated snakemake command line.
>
> Default (same target as full workflow): `snakemake --report results/report.html`.
>
> Specify targets (output files or rules, separated by a space): `snakemake multiqc parse_fasta --report results/report.html`

> [!TIP]
> 
> If one were to generate their own (perfectly formatted) file, for example the list of expected mutants to accomodate a custom experimental design, it is possible to plug it using `--touch`
>
> One compressed CSV file per mutated locus is expected in the `project_files/expected_mut` folder (create it if it does not exist)
> and should be named using the same values as in [the dataframe containing wild-type sequences](../config/project_files/wt_seq.csv), 
> e.g. `project_files/expected_mut/CN_F1.csv.gz`
>
> The following: `snakemake --touch config/project_files/expected_mut/CN_F1.csv.gz` instructs the workflow not to overwrite the file provided by the user and use it as if it had been generated by the workflow

## Apptainer support

> [!NOTE]
> 
> Apptainer is currently not supported... although it might be in the future!

Run the workflow using: `snakemake --profile profiles/slurm --sdm conda apptainer`. The container should be created first, then conda envs will be created for each rule inside the container. This option is meant to be used on a system where you want to isolate the (many) files installed by `conda`. This option is **not** suited for local execution.

## Edit pipeline

> [!IMPORTANT]
> 
> On can only modify the pipeline after cloning the repo, not upon snakedeployment.


One can manually edit the [Snakefile](../workflow/Snakefile) and/or the rules (.smk files in rules folder) to edit the main steps of the pipeline. This should not be required to run the standard pipeline and should be done only when the core workflow itself needs to be modified.

> [!TIP]
> 
> In certain cases, it might be interesting to modify the scripts themselves, for example one might want to alter **plotting**.
> The recommended way is to edit [the functions](../workflow/scripts/plotting_functions.py) and call them from a custom script.

## Edit Jupyter notebooks

What follows is related to the fact that earlier versions of gyōza relied on Jupyter notebooks instead of Python scripts. I'm keeping this here if anyone working on a Snakemake workflow is interested.
    
Manual retrieval and modification of template notebooks is **not recommended** (you would need to activate the proper conda env or any env with the same packages installed, then export the kernelspecs, then modify the paths and parameters by hand, which can lead to some typing errors, etc)

Thankfully, there is a snakemake command that allows **interactive editing of any template notebook**, using any output file (from the notebook) as argument. The following example will generate URLs to open `jupyter`, in which we can edit the process_read_counts notebook that outputs the upset_plot.svg file (this example is unfortunately no longer valid, since this specific output file is now generated by a Python script)

```
snakemake --use-conda --cores 1 --edit-notebook results/graphs/upset_plot.svg
```

> [!WARNING]
> 
> If you are running `snakemake` on a server, you might need to open a SSH tunnel between your local machine and the server by running the following command from a local terminal (should not be necessary when running locally on your machine): `ssh -L 8888:localhost:8888 <USER>@<ADRESS>`
>
> Adapt port if necessary, 8888 or 8889, should match what is featured in the URL generated by snakemake/jupyter

The command simply opens a terminal on the server, but now you can copy-paste the URL in your local browser.
    
You can then open the notebook, run it (kernel and paths taken care of) and save it once the modifications have been done. Then click on "Close notebook" and finally "Shut down" on the main jupyter dashboard. The quitting and saving should be displayed on the initial server's terminal (the one from which the `snakemake` command was run). You'll also notice that the size of the template notebook file could be smaller, because outputs are automatically flushed out. To retrieve a notebook with outputs (for future HTML export for example), locate the notebook in the appropriate folder once the pipeline has run (path specified in the log directive of the .smk file).