Installation
============

MD-SAPT can be installed from the psi4 Conda repo like so:

.. code-block:: bash

    conda install -c psi4/label/dev -c conda-forge mdsapt

Alternatively, it can be installed by cloning the GitHub repository.

.. code-block:: bash

    git clone https://github.com/calpolyccg/MDSAPT.git
    pip install ./MDSAPT

To ensure it's been installed correctly, run `mdsapt` or `python3 -m mdsapt`.

.. code-block::
    Warning: importing 'simtk.openmm' is deprecated.  Import 'openmm' instead.
    2022-03-30 09:32:50,071 mdsapt       INFO     MDSAPT 1.2.0 starting
    2022-03-30 09:32:50,071 mdsapt       INFO     Copyright (c) 2021 Alia Lescoulie, Astrid Yu, and Ashley Ringer McDonald
    2022-03-30 09:32:50,071 mdsapt       INFO     Released under GPLv3 License
    Usage: python -m mdsapt [OPTIONS] COMMAND [ARGS]...

      MDSAPT - Molecular Dynamics Symmetry-Adapted Perturbation Theory

      This command-line interface lets you easily do common MDSAPT-related tasks.

    Options:
      --help  Show this message and exit.

    Commands:
      generate  Generate a template input file at filename.
      run       Run a SAPT calculation using the configuration in in_file.
