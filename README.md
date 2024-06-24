MD-SAPT
==============================
SAPT Calculations for MDAnalysis

[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/calpolyccg/MDSAPT/workflows/CI/badge.svg)](https://github.com/alescolie/MDSAPT/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/calpolyccg/MDSAPT/branch/master/graph/badge.svg)](https://codecov.io/gh/alescoulie/MDSAPT/branch/master)
[![Documentation Status](https://readthedocs.org/projects/mdsapt/badge/?version=latest)](https://mdsapt.readthedocs.io/en/latest/?badge=latest)
[![Launch binder demo](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/calpolyccg/MDSAPT_demo/master?labpath=MD-SAPT_demo.ipynb)

An MDAnalysis-kit for calculating SAPT of molecular dynamics trajectories in psi4. [Click here for a demo!](https://mybinder.org/v2/gh/calpolyccg/MDSAPT_demo/master?labpath=MD-SAPT_demo.ipynb)

### Installation

Installing with Conda on Python version 3.9 to 3.11 on either MacOS (arm or x86) or Linux using the command:

``` bash
conda install -c conda-forge mdsapt
```

MD-SAPT can also be installed from source:

``` bash
git clone https://github.com/calpolyccg/MDSAPT.git
pip install ./MDSAPT
```

If you use [nix](https://nixos.wiki/) MD-SAPT can be also be build with:

```
git clone https://github.com/calpolyccg/MDSAPT.git
cd MDSAPT
nix build
```

You can also incorporate it to a project managed with nix by adding it to the imports of your `flake.nix` or `shell.nix`.

### Contributing
A development environment can be created using nix:

``` bash 
git clone https://github.com/calpolyccg/MDSAPT.git

cd MDSAPT

nix develop
```

### Copyright

Copyright (c) 2021-2024, Alia Lescoulie

#### Acknowledgments

This work was supported by the Bill and Linda Frost Fund at Cal Poly San Luis Obispo.

We acknowledge the support of NSF award CHE-2018427.  This award provided the computational resources to support this project through the MERCURY Consortium Skylight cluster on Palmetto. 

Project based on the [Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.

