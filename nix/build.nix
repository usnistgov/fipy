{ nixpkgs, pypkgs, pysparse, preshellhook }:
let
  gmsh = import ./gmsh.nix { inherit nixpkgs; };
  skfmm = import ./skfmm.nix { inherit pypkgs; };
in
  pypkgs.buildPythonPackage rec {
     pname = "fipy";
     version = "3.1.3.dev";
     env = nixpkgs.buildEnv { name=pname; paths=buildInputs; };
     buildInputs = [
       pypkgs.python
       pypkgs.numpy
       pypkgs.scipy
       gmsh
       skfmm
       pysparse
       pypkgs.matplotlib
       pypkgs.tkinter
     ];
     src=./..;
     doCheck=false;
     preShellHook = preshellhook version;
     prePatch = preshellhook version;
     meta = {
       homepage = "https://www.ctcms.nist.gov/fipy/";
       description = "A Finite Volume PDE Solver Using Python";
       version = version;
       license = nixpkgs.stdenv.lib.licenses.free;
     };
     catchConflicts=false;
  }
