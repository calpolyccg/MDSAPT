Quickstart Guide
================

MD-SAPT simplifies the calculation of SAPT interaction energies between selected residue pairs in MD trajectories. Running it just requires MD simulation files.

Prerequisites
_____________

Ensure that you have the following things set up:

 - You have existing MD trajectory and topology files in any form that `MDAnalysis <https://mdanalysis.org>`_ You have
already installed MDSAPT, following the :ref:`installation page <install>`.

.. note:
    If your `PATH` environment variable is not set up to point to installed Python modules, then invoking `mdsapt` directly, as shown in this guide, may not work. In that case, try running `python3 -m mdsapt` instead.

.. _gen-input:
Generating an input file
________________________

The following steps describe how to set up the input YAML file.

 1. Run ``mdsapt generate [filename]`` to generate a blank Trajectory SAPT input file at the given filename. Don't
worry, it will not overwrite any files unless you explicitly provide the `-f` flag.

 2. Specify the MD file locations in the input file

 3. Specify the residue numbers you wish to analyze

 4. Specify the pairs of residues from the numbers in step 3

 5. Specify trajectory, optimization, and SAPT settings

Here is an example of a filled-out YAML file:

.. code-block:: yaml

    psi4:
      method: "sapt0"
      basis: "jun-cc-pvdz"
      settings:
        reference: "rhf"
      save_output: true
    simulation:
      ph: 7.0
      charge_guesser: "standard"
      # charge_guesser: 'rdkit'  # to use rdkit. Make sure it is installed first.
    system_limits:
      ncpus: 32
      memory: "80GB"
    analysis:
      ### This section is for running TrajectorySAPT. To run other types of analyses, see below.
      type: "trajectory"

      topology: testtop.psf
      trajectories:
        - testtraj.dcd
      pairs:
        # Place pair of  selections defined above in a list of lists
        - [109, 196]
        - [197, 199]
        - [208, 200]
        - [156, 44]
        - [84, 13]
      frames:
        start: 78
        stop: 97
        step: 1
      output: "output.csv"


Running SAPT
____________

Using the `mdsapt` CLI
^^^^^^^^^^^^^^^^^^^^^^

This mode is suited for:

 - running on an HPC cluster
 - running a simple input file

With the input done MD-SAPT is ready to be run. Simply execute

.. code-block:: bash
   mdsapt run [filename] [output]

and it will run SAPT on your trajectory using the parameters specified in your input file.

Using the `mdsapt` Python library
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This mode is suited for:

 - using MD-SAPT in a notebook
 - using MD-SAPT in your own library or applications

The classes involved are as follows:

 - The settings are read using :class:`mdsapt.reader.InputReader`
 - The `InputReader` is then passed into :class:`mdsapt.reader.Optimizer` which handles preparing residues.
 - Finally, :class:`mdsapt.reader.TrajectorySAPT` is used to run SAPT over the MD data.
 - The results are stored in a :class:`Pandas.DataFrame` which can be accessed under the `TrajectorySAPT.results` property.

Here is some code demonstrating it:

.. code-block:: Python

    import mdsapt


    config = mdsapt.load_from_yaml_file('runinput.yaml')
    sapt_run = mdsapt.TrajectorySAPT(config)
    sapt_run.run(config.start, config.stop, config.step)
    sapt_run.results.to_csv('results.csv')

See also `the Binder demo <https://mybinder.org/v2/gh/calpolyccg/MDSAPT_demo/master?labpath=MD-SAPT_demo.ipynb>`_ for a bigger example.

