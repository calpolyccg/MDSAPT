Quickstart Guide
================

MD-SAPT simplifies for the calculation of SAPT interaction energies from MD
trajectories. Running it just requires MD simulation files.

Basic work flow
===============

MD-SAPT can be run as a script or in a notebook. It gives the SAPT energy
between the specified residue pairs.

Setup Process
_____________

The following steps describe how to set up the input yaml file

 1. Use the `mdsapt_get_runinput.py` script to generate a blank input file

 2. Specify the MD file locations in the input file

 3. Specify the residue numbers you wish to analyze

 4. Specify the pairs of residues from the numbers in step 3

 5. Specify trajectory, optimization, and SAPT settings

Running SAPT
____________

With the input done MD-SAPT is ready to be run. The settings are read using :class:`mdsapt.reader.InputReader`  which is then passed into :class:`mdsapt.reader.Optimizer` which handles preparing residues. Finally :class:`mdsapt.reader.TrajectorySAPT` is used to run SAPT over the MD data. The results are stored in a :class:`Pandas.DataFrame`.

.. code-block:: Python

    import mdsapt

    Settings = mdsapt.InputReader('runinput.yaml')

    Opt = mdsapt.Optimizer(Settings)

    SAPT_run = mdsapt.TrajectorySAPT(Settings, Opt)

    SAPT_run.run(Settings.start, Settings.stop, Settings.step)

    SAPT_run.results.to_csv('results.csv')