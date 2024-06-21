{
  description = "SAPT energy calculator built using MDAnalysis and Psi4 ";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
    poetry2nix.url = "github:nix-community/poetry2nix";
    qchem.url = "github:Nix-QChem/NixOS-QChem";
  };

  outputs = { self, nixpkgs, flake-utils, poetry2nix, qchem }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        lib = nixpkgs.lib;

        pkgs = import nixpkgs {
          system = system;
          overlays = [ 
            qchem.overlays.qchem 
            poetry2nix.overlays.default
          ];
        };


        # Build dependencies for pacakges
        pypkgs-build-reqs = {
          mdanalysis = [ "setuptools" ];
          mda-xdrlib = [ "setuptools" ];
          mmtf-python = [ "setuptools" ];
          mrcfile = [ "setuptools" ];
          griddataformats = [ "setuptools" ];
          pyright = [ "setuptools" ];
          #matplotlib = [ "pybind11" ];
          #scipy = [ "setuptools" "wheel" "pybind11" "pythran" ];
        };

        p2n-overrides = pkgs.poetry2nix.defaultPoetryOverrides.extend (final: prev:
          builtins.mapAttrs (package: build-reqs:
            (builtins.getAttr package prev).overridePythonAttrs (old: {
              buildInputs = (old.buildInputs or [ ]) ++ (builtins.map (pkg: if builtins.isString pkg then builtins.getAttr pkg prev else pkg) build-reqs);
            })
          ) pypkgs-build-reqs
        );

        unfuckScipy = final: prev: {
          scipy = pkgs.python311Packages.scipy;
        };

        unfuckNumpy = final: prev: {
          numpy = pkgs.python311Packages.numpy;
        };

        python = pkgs.python311.override {
          packageOverrides = unfuckScipy;
        };
        
        poetryEnv = pkgs.poetry2nix.mkPoetryEnv {
          projectDir = ./.;
          python = python;
          overrides = [ p2n-overrides unfuckScipy unfuckNumpy ]; #unfuckScipy ];
        };

        devEnv = pkgs.mkShell {
          propagatedBuildInputs = [ poetryEnv];
          buildInputs = with pkgs; [pkgs.qchem.psi4 pkgs.qchem.openmm pkgs.qchem.pdbfixer pkgs.act ];
        };

        # DON'T FORGET TO PUT YOUR PACKAGE NAME HERE, REMOVING `throw`
        packageName = "MD-SAPT";

      in {
        devShells.default = devEnv;
        
        test-ci = pkgs.runCommand "test-ci" {} ''
          act
        '';

        test-mdsapt = pkgs.runCommand "test-mdsapt" {} ''
          pytest -v ./mdsapt --cov=mdsapt --cov-report=xml
        '';

      });
}

