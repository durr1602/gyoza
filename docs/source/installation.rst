Installation
============

.. tabs::

   .. tab:: Linux

      Open a terminal.

   .. tab:: Mac OSX

      Open a terminal.

   .. tab:: Windows

      Install WSL2 (`WSL installation instructions <https://learn.microsoft.com/en-us/windows/wsl/install>`__)
      , then open a WSL2 terminal.

1. Make sure you have ``conda>=24.9.1`` (see `Conda documentation <https://conda.io/docs/index.html>`__).
Proceed to step 2 or 3 accordingly.
2. To update conda: ``conda update -n base -c defaults conda --repodata-fn=repodata.json``, or ``mamba update conda``.
If you have multiple instances of conda and want to uninstall, say, Anaconda, run ``rm -rf anaconda3``.
3. To install Miniforge (conda), follow the instructions listed in the “Step **1b**” section of `this
tutorial <https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-1b-installing-miniforge>`__.
4. Install ``uv``: `https://github.com/astral-sh/uv <https://github.com/astral-sh/uv>`__
5. Create a Python>=3.12 virtual environment to install dependencies.
The specified Python version should automatically be installed if not available locally.

::

   uv venv gyoza_env --python 3.13
   source gyoza_env/bin/activate
   uv pip install "snakedeploy>=0.11.0" "snakemake>=9.9.0" "snakemake-wrapper-utils>=0.7.2" pygments snakemake-executor-plugin-cluster-generic setuptools

Whenever your gyoza_env is activated, you should see it in the prompt:

::

   (gyoza_env) <USER>@<MACHINE>

If another env is displayed, e.g. (base), make sure to deactivate it so
that there is only (gyoza_env) left. You may need to run
``conda config --set auto_activate_base False``.

Recommended installation of gyōza
---------------------------------

Deploy gyōza
~~~~~~~~~~~~

The following command line uses
`Snakedeploy <https://snakedeploy.readthedocs.io/en/latest/index.html>`__
to create a minimal file tree that will allow you to use gyōza.

It is specific to both the version of gyōza deployed and the DMS project
your want to analyze, and should make it easier to create a repository
for improved reproducibility. The ``--tag`` argument accepts any branch or
release version tag

::

   snakedeploy deploy-workflow https://github.com/durr1602/gyoza my_gyoza_project --tag main
   cd my_gyoza_project

Alternative option
------------------

Clone this repository
~~~~~~~~~~~~~~~~~~~~~

This option provides the most flexibility, although should be reserved for advanced users.

::

   git clone https://github.com/durr1602/gyoza.git
   cd gyoza
