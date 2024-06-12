{
  description = "SAPT energy calculator built using MDAnalysis and Psi4 ";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-23.11";
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

      
        cython = with pkgs.python3Packages; cython_3.overrideAttrs (oldAttrs: rec {
          version = "3.0.10";
          src = fetchPypi {
            pname = "Cython";
            inherit version;
            hash = "sha256-3MlnOTMfuFTc9QP5RgdXbP6EiAZsYcpQ39VYNvEy3pk=";
          };
        });

        pypkgs-build-reqs = {
          mdanalysis = [ "setuptools" ];
          mda-xdrlib = [ "setuptools" ];
          mmtf-python = [ "setuptools" ];
          mrcfile = [ "setuptools" ];
          scipy = [ "setuptools" "wheel" "pybind11" "pythran" ];
        };

        p2n-overrides = pkgs.poetry2nix.defaultPoetryOverrides.extend (final: prev:
          builtins.mapAttrs (package: build-reqs:
            (builtins.getAttr package prev).overridePythonAttrs (old: {
              buildInputs = (old.buildInputs or [ ]) ++ (builtins.map (pkg: if builtins.isString pkg then builtins.getAttr pkg prev else pkg) build-reqs);
            })
          ) pypkgs-build-reqs
        );

        unfuckScipy = final: prev: {
          cython = cython;
          cython_3 = cython;
          scipy = prev.scipy.overridePythonAttrs (old: {
            buildInputs = old.buildInputs ++(with final; [ cython setuptools wheel pythran pybind11 ]);
          });
        };

        devEnv = pkgs.poetry2nix.mkPoetryEnv {
          projectDir = ./.;
          extraPackages = (ps: [ pkgs.qchem.psi4 pkgs.qchem.openmm ]);
          python = pkgs.python311.override {
            packageOverrides = pfinal: pprev: {
              cython = cython;
              cython_3 = cython;
            };
          };
          overrides = [p2n-overrides unfuckScipy ];
        };

        # DON'T FORGET TO PUT YOUR PACKAGE NAME HERE, REMOVING `throw`
        packageName = "MD-SAPT";

      in {
        test.cython = cython;
        devShells.default = devEnv;
      });
}

