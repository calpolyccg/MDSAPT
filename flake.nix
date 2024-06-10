{
  description = "SAPT energy calculator built using MDAnalysis and Psi4 ";

  inputs = {
    nixpkgs.url = "github:nixos/nixpkgs/nixos-unstable";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils}:
    flake-utils.lib.eachDefaultSystem (system:
      let
        lib = nixpkgs.lib;
        pkgs = nixpkgs.legacyPackages.${system};



        # DON'T FORGET TO PUT YOUR PACKAGE NAME HERE, REMOVING `throw`
        packageName = "MD-SAPT";

      in {

        devShells.default = pkgs.mkShell {
          buildInputs = [ 
            pkgs.poetry
          ];
        };
      });
}

