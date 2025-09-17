Installation
============

If you use Windows, the following steps need to be run in WSL2 (`WSL
installation
instructions <https://learn.microsoft.com/en-us/windows/wsl/install>`__).

If you already have WSL2 or if you use Linux or MacOS X, the following
steps need ``conda`` (see `Conda
documentation <https://conda.io/docs/index.html>`__).

If you don’t already have it, we recommend installing Miniforge by
following the instructions listed in the “Step **1b**” section of `this
tutorial <https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-1b-installing-miniforge>`__.

If you already have it, make sure to update it to a recent version
(>=24.7.1, ideally even more recent such as >=24.9.1). Here are two
helpful command lines to do this:
``conda update -n base -c defaults conda --repodata-fn=repodata.json``,
or ``mamba update conda``. If you have multiple instances of conda and
want to uninstall, say, Anaconda, run ``rm -rf anaconda3``.

1. Install ```uv`` <https://github.com/astral-sh/uv>`__
2. Create a Python>=3.12 virtual environment to install dependencies.
   The specified Python version should automatically be installed if not
   available locally.

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

2. Deploy gyōza.

The following command line uses
`Snakedeploy <https://snakedeploy.readthedocs.io/en/latest/index.html>`__
to create a minimal file tree that will allow you to use gyōza.

It is specific to both the version of gyōza deployed and the DMS project
your want to analyze, and should make it easier to create a repository
for improved reproducibility. The –tag argument accepts any branch or
release version tag

::

   snakedeploy deploy-workflow https://github.com/durr1602/gyoza my_gyoza_project --tag main
   cd my_gyoza_project

Alternative option: clone the repository
----------------------------------------

This option provides the most flexibility, since the entire repository
is cloned, meaning one can potentially modify the code of gyōza.

2. Clone this repository:

::

   git clone https://github.com/durr1602/gyoza.git
   cd gyoza
