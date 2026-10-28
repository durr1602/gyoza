Installation
============

Open a terminal
---------------

.. tabs::

    .. tab:: Linux

       Open a terminal.

    .. tab:: Mac OSX

       If your Mac is equipped with an Intel processor, simply open a terminal.

       If your Mac is equipped with Apple silicon (M1/M2/M3), you may need to
       `install Rosetta <https://support.apple.com/en-ca/102527>`__ first, then
       open a terminal.

    .. tab:: Windows

       Install WSL2 (`WSL installation instructions <https://learn.microsoft.com/en-us/windows/wsl/install>`__)
       , then open a WSL2 terminal.

Install requirements
--------------------

Snakemake relies on ``conda`` to create the software environment of each step. Choose
one of the following options.

.. tabs::

    .. tab:: pixi (no conda install)

       `pixi <https://pixi.sh>`__ installs ``conda`` inside the ``gyoza_env``
       environment only, nothing is installed system-wide.

       1. Install ``pixi``: https://pixi.sh/latest/installation/
       2. Create the ``gyoza_env`` environment:

       ::

           pixi init gyoza_env
           pixi add --manifest-path gyoza_env "python=3.13" "conda>=24.9.1"
           pixi add --manifest-path gyoza_env --pypi "snakedeploy==0.15.0" "snakemake==9.14.0" "snakemake-wrapper-utils==0.8.0" "snakemake-executor-plugin-cluster-generic==1.0.9" pygments setuptools

       3. Activate it (``exit`` to deactivate):

       ::

           pixi shell --manifest-path gyoza_env

    .. tab:: conda + uv

       1. Make sure you have ``conda>=24.9.1`` (see `Conda documentation
          <https://conda.io/docs/index.html>`__). Proceed to step 2 or 3 accordingly.
       2. To update ``conda``: ``conda update -n base -c defaults conda
          --repodata-fn=repodata.json``, or ``mamba update conda``. If you have multiple
          instances of ``conda`` and want to uninstall, say, Anaconda, run ``rm -rf anaconda3``.
       3. To install Miniforge (``conda``), follow the instructions listed in the “Step **1b**\ ”
          section of `this tutorial
          <https://snakemake.readthedocs.io/en/stable/tutorial/setup.html#step-1b-installing-miniforge>`__.
       4. Install ``uv``: https://github.com/astral-sh/uv
       5. Create a Python>=3.12 virtual environment to install dependencies. The specified
          Python version should automatically be installed if not available locally.

       ::

           uv venv gyoza_env --python 3.13
           source gyoza_env/bin/activate
           uv pip install "snakedeploy==0.15.0" "snakemake==9.14.0" "snakemake-wrapper-utils==0.8.0" "snakemake-executor-plugin-cluster-generic==1.0.9" pygments setuptools

       If another env is displayed in the prompt, e.g. ``(base)``, make sure to
       deactivate it so that there is only ``(gyoza_env)`` left. You may need to run
       ``conda config --set auto_activate_base False``.

Whenever ``gyoza_env`` is activated, you should see it in the prompt:

::

    (gyoza_env) <USER>@<MACHINE>

Get gyōza
---------

.. tabs::

    .. tab:: Snakedeploy (recommended)

       The following command line uses `Snakedeploy
       <https://snakedeploy.readthedocs.io/en/latest/index.html>`__ to create a minimal file
       tree that will allow you to use gyōza.

       It is specific to both the version of gyōza deployed and the DMS project your want to
       analyze, and should make it easier to create a repository for improved reproducibility.
       Specify a release version tag for the ``--tag`` argument (latest release:
       |latest_tag|), or a specific commit to use a development version (advanced).

       .. parsed-literal::

           snakedeploy deploy-workflow https://github.com/durr1602/gyoza my_gyoza_project --tag |latest_tag|
           cd my_gyoza_project

    .. tab:: Clone the repository (advanced)

       This option provides the most flexibility, although should be reserved for advanced
       users.

       ::

           git clone https://github.com/durr1602/gyoza.git
           cd gyoza
