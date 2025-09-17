Usage
=====

Prepare files and edit config
-----------------------------

.. important::
   1. Read the `config documentation <../config/README.md>`__
   2. **Edit** `the main config <../config/config.yaml>`__

   If you choose to execute gyōza locally, you may optionally edit `the
   default profile <../profiles/slurm/config.v8+.yaml>`__

   If you choose to execute gyōza with SLURM, don’t forget to edit `the
   slurm profile <../profiles/slurm/config.v8+.yaml>`__

(optional) Prepare environments
-------------------------------

.. tip::
   3. (optional) Create all conda environments using:
      ``snakemake --conda-create-envs-only``.

   This step is only required before first use and is always included
   when running the workflow. Unfortunately, at the time of writing,
   validations will be run even for these command lines. This means that
   you need to fully prepare the workflow before creating all envs.

Check pipeline
--------------

.. important::
   4. Perform a dry run using: ``snakemake -n``

   This step is strongly recommended. It will make sure the prepared
   workflow does not contain any error and will display the rules
   (steps) that need to be run in order to run the workflow (built
   dynamically based on the config). If you’re running the workflow for
   the first time and you toggled in normalization with growth data, you
   should see a warning prompting you to edit the generated template
   file (for more details, go back to step 1).

Run pipeline
------------

.. important::
   5. Run the workflow either **locally (a)** or **using SLURM (b)**

   a) Local execution: ``snakemake`` or
      ``snakemake --profile profiles/default``
   b) SLURM (1 job per rule per sample):
      ``snakemake --profile profiles/slurm``. Jobs wait in the queue
      until the resources are allocated. For example, if you’re allowed
      40 CPUs, only 4 jobs at 10 CPUs each will be able to run at once.
      Once those jobs are completed, the next ones in the queue will
      automatically start. Fore more info on cluster execution: read the
      doc on `smk-cluster-generic
      plugin <https://github.com/jdblischak/smk-simple-slurm/tree/main>`__

Abort pipeline / exit terminal
------------------------------

If snakemake is launched directly from the command line, the process
will be output to the terminal. Exiting with ``<Ctrl+C>`` is currently
interpreted (as specified in the `the slurm
profile <../profiles/slurm/config.v8+.yaml>`__) as cancelling all
submitted jobs (``scancel``). Exiting during a local execution will
**also** abort the workflow. This means that while the workflow is
running, the user cannot get the prompt back.

There are 3 possible options to get the prompt back and/or exit the
terminal without aborting the workflow: 1. Open a new tab on your
terminal (may require to log into the session again) 2. Use ``nohup``
(e.g. ``nohup snakemake``). Closing the tab will not abort the workflow.
3. Use the terminal multiplexer
```tmux`` <https://github.com/tmux/tmux/wiki/Getting-Started>`__. This
way you can get the prompt back right away, exit without aborting and
even reconnect to the session from a different machine.

Please make sure ``tmux`` is installed (already installed on some
servers). Then, follow the steps: 1. Type ``tmux new -s snakes`` to
launch a new tmux session 2. Make sure the gyoza_env is activated (or
activate it) 3. Make sure you are in the proper directory and launch the
pipeline with ``snakemake --profile profiles/slurm`` 4. To close
(detach) the session, type ``<Ctrl+b>``, then ``<d>``. You should see
the message: ``[detached (from session snakes)]`` 5. To reconnect
(attach) to the session, for example from a different machine:
``tmux attach -t snakes``. You can also see existing sessions with
``tmux ls``. 6. To close the session when everything is finished, type
``<Ctrl+b>``, then ``<:>``, then ``kill-session`` and finally
``<Enter>``.

Useful Snakemake command lines
------------------------------

.. tip::
   For both the dry run and the actual run, you can decide to run the
   workflow only until a certain file is generated or rule is completed,
   using the ``--until`` flag in the snakemake command line

   For example: ``snakemake -n --until stats``

.. tip::
   Similarly, you can omit rules (and depending downstream rules) from
   the workflow.

   For example: ``snakemake --omit-from rule_to_omit``

.. note::
   To account for intermediate runs, the report is built dynamically
   based on existing files (out of all files expected to be included in
   the report). To manually specify what should be included in the
   report, use the dedicated snakemake command line.

   Default (same target as full workflow):
   ``snakemake --report results/report.html --report-stylesheet config/style/report-stylesheet.css``.

   Specify targets (output files or rules, separated by a space):
   ``snakemake multiqc parse_fasta --report results/report.html --report-stylesheet config/style/report-stylesheet.css``

.. tip::
   *This particular tip might be obsolete given recent improvements, but
   I’m leaving it here for now*

   If one were to generate their own (perfectly formatted) file, for
   example the list of expected mutants to accomodate a custom
   experimental design, it is possible to plug it using ``--touch``

   One compressed CSV file per mutated locus is expected in the
   ``project_files/expected_mut`` folder (create it if it does not
   exist) and should be named using the same values as in `the dataframe
   containing wild-type
   sequences <../config/project_files/wt_seq.csv>`__,
   e.g. ``project_files/expected_mut/CN_F1.csv.gz``

   The following:
   ``snakemake --touch config/project_files/expected_mut/CN_F1.csv.gz``
   instructs the workflow not to overwrite the file provided by the user
   and use it as if it had been generated by the workflow

Apptainer support
-----------------

.. note::
   Apptainer is currently not supported… although it might be in the
   future!

Run the workflow using:
``snakemake --profile profiles/slurm --sdm conda apptainer``. The
container should be created first, then conda envs will be created for
each rule inside the container. This option is meant to be used on a
system where you want to isolate the (many) files installed by
``conda``. This option is **not** suited for local execution.

Edit pipeline
-------------

.. important::
   On can only modify the pipeline after cloning the repo, not upon
   snakedeployment.

One can manually edit the `Snakefile <../workflow/Snakefile>`__ and/or
the rules (.smk files in rules folder) to edit the main steps of the
pipeline. This should not be required to run the standard pipeline and
should be done only when the core workflow itself needs to be modified.

.. tip::
   In certain cases, it might be interesting to modify the scripts
   themselves, for example one might want to alter **plotting**. The
   recommended way is to edit `the
   functions <../workflow/scripts/plot_scores.py>`__ and call them from
   a custom script.
