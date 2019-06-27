# Usage:
#
#     $ nix-shell --pure --argstr python_version 36
#
# to use Python version 3.6 for example. Tested with 2.7, 3.6 and 3.7.


{ nixpkgs ? (import ./nix/nixpkgs_version.nix), python_version ? "27" }:
let
  name = "python" + python_version + "Packages";
  pypkgs = nixpkgs.${name};
  pysparse = if (python_version == "27") then
    (import ./nix/pysparse.nix { inherit nixpkgs pypkgs; })
  else
    null;
  gmsh = import ./nix/gmsh.nix { inherit nixpkgs; };
  skfmm = import ./nix/skfmm.nix { inherit nixpkgs pypkgs; };
in
  pypkgs.buildPythonPackage rec {
    name = "fipy";
    env = nixpkgs.buildEnv { name=name; paths=buildInputs; };
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
      pypkgs.future
    ];
    src=builtins.filterSource (path: type: type != "directory" || baseNameOf path != ".git") ./.;
    doCheck=false;
    meta = {
      homepage = "https://www.ctcms.nist.gov/fipy/";
      description = "A Finite Volume PDE Solver Using Python";
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
