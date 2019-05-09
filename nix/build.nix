{ nixpkgs, pypkgs, pysparse, preshellhook }:
let
  gmsh = import ./gmsh.nix { inherit nixpkgs; };
  skfmm = import ./skfmm.nix { inherit nixpkgs pypkgs; };
in
  pypkgs.buildPythonPackage rec {
     pname = "fipy";
     version = "3.1.3.dev";
     env = nixpkgs.buildEnv { name=pname; paths=buildInputs; };
     buildInputs = [
       pypkgs.pip
       pypkgs.python
       pypkgs.numpy
       pypkgs.scipy
       gmsh
       skfmm
       pysparse
       pypkgs.matplotlib
       pypkgs.tkinter
       nixpkgs.pkgs.git
       nixpkgs.imagemagick
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
     postShellHook = ''
       SOURCE_DATE_EPOCH=$(date +%s)
       export PYTHONUSERBASE=$PWD/.local
       export USER_SITE=`python -c "import site; print(site.USER_SITE)"`
       export PYTHONPATH=$PYTHONPATH:$USER_SITE
       export PATH=$PATH:$PYTHONUSERBASE/bin

       ## To build the docs
       # pip install --user sphinx
       # pip install --user sphinxcontrib-bibtex
       # pip install --user git+https://github.com/thewtex/sphinx-contrib.git#subdirectory=traclinks
       # required for embedded plots in documentation
       # pip install --user pandas
       # python setup.py build_docs --html

       ## To build PyAMG add nixpkgs.gcc to buildInputs and then
       # pip install --user pyamg
     '';
  }
