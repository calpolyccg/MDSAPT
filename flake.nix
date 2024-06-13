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

        # Version overrides for scipy 
        /*
        cython = with pkgs.python311Packages; cython.overrideAttrs (oldAttrs: rec {
          version = "3.0.10";
          src = fetchPypi {
            pname = "Cython";
            inherit version;
            hash = "sha256-3MlnOTMfuFTc9QP5RgdXbP6EiAZsYcpQ39VYNvEy3pk=";
          };
        });
        */



        #pybind11-212 = with pkgs.python311Packages; pybind11.overrideAttrs (oldAttrs: rec {
        #  version = "2.12.0";
        #  src = fetchPypi {
        #    pname = "pybind11";
        #    inherit version;
        #    hash = "sha256-XjxVeoSwa5aSR2MEB/xNmFvtFXtCU7ExU7jo4WXgw9w=";
        #  };
        #});

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
          /*
          scipy = prev.scipy.overridePythonAttrs (old: {
            nativeBuildInputs = with final; [cython_3];
            buildInputs = old.buildInputs ++ (with final; [ cython_3 setuptools wheel pythran pybind11 ]);
          });
          */
          scipy = pkgs.python311Packages.scipy;
        };

        unfuckNumpy = final: prev: {
          numpy = pkgs.python311Packages.numpy;
        };

        python = pkgs.python311.override {
          packageOverrides = unfuckScipy;
        };
        
        #unfuckMatPlotLib = final: prev: {
        #  pybind11 = pybind11-212;
        #  matplotlib = prev.matplotlib.overridePythonAttrs (old: {
        #    buildInputs = old.buildInputs ++(with final; [ pybind11 ]);
        #  });
        #};

        poetryEnv = pkgs.poetry2nix.mkPoetryEnv {
          projectDir = ./.;
          python = python;
          overrides = [ p2n-overrides unfuckScipy unfuckNumpy ]; #unfuckScipy ];
        };

        devEnv = pkgs.mkShell {
          propagatedBuildInputs = [ poetryEnv];
          buildInputs = with pkgs; [pkgs.qchem.psi4 pkgs.qchem.openmm pkgs.qchem.pdbfixer ];
        };

        # DON'T FORGET TO PUT YOUR PACKAGE NAME HERE, REMOVING `throw`
        packageName = "MD-SAPT";

      in {
        devShells.default = devEnv;
      });
}

