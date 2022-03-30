Quickstart Guide
================

MD-SAPT simplifies the calculation of SAPT interaction energies between selected residue pairs in MD trajectories. Running it just requires MD simulation files.

Prerequisites
_____________

Ensure that you have the following things set up:

 - You have existing MD trajectory and topology files in any form that `MDAnalysis <https://mdanalysis.readthedocs.io/en/latest/`_ supports
 - You have already installed MDSAPT, following the :installation guide:`install`

.. note:
    If your `PATH` environment variable is not set up to point to installed Python modules, then invoking `mdsapt` directly, as shown in this guide, may not work. In that case, try running `python3 -m mdsapt` instead.

Generating an input file
________________________

The following steps describe how to set up the input YAML file.

 1. Run `mdsapt generate [filename]` to generate a blank Trajectory SAPT input file at the given filename. Don't worry, it will not overwrite any files unless you explicitly provide the `-f` flag.

 2. Specify the MD file locations in the input file

 3. Specify the residue numbers you wish to analyze

 4. Specify the pairs of residues from the numbers in step 3

 5. Specify trajectory, optimization, and SAPT settings

Here is an example of a filled-out YAML file:

.. code-block:: yaml
    topology_path: 'testtop.psf'
    trajectory_paths:
      - 'testtraj.dcd'
    selection_resid_num:
      - 11
      - 199
    int_pairs:
      # Place pair of  selections defined above in a list of lists
      - [11, 199]
    trajectory_settings:
      start: 0
      stop: 98
      step: 1
    system_settings:
      ncpus: '12'
      memory: '12GB'
      time: '24:00:00'
    opt_settings:
      pH: 7
    sapt_settings:
      method: 'sapt0'
      basis: 'jun-cc-pvdz'
      settings:
        reference: 'rhf'
      save_psi4_output: true

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


    Settings = mdsapt.InputReader('runinput.yaml')
    Opt = mdsapt.Optimizer(Settings)
    SAPT_run = mdsapt.TrajectorySAPT(Settings, Opt)
    SAPT_run.run(Settings.start, Settings.stop, Settings.step)
    SAPT_run.results.to_csv('results.csv')

See also `the Binder demo <https://mybinder.org/v2/gh/calpolyccg/MDSAPT_demo/master?labpath=MD-SAPT_demo.ipynb>`_ for a bigger example.

