{ nixpkgs ? import ./nix/nixpkgs_version.nix }:
let
  pypkgs = nixpkgs.python36Packages;
  pysparse = null;
  preshellhook = version: ''
    echo hello
    2to3 --write . &> /dev/null;
    2to3 --write --doctests_only . &> /dev/null;
    sed -i 's/version = getVersion()/version="${version}"/' setup.py
  '';
in
  import ./build.nix { inherit nixpkgs pypkgs pysparse preshellhook; }

## After running nix-shell
##
## $ find . -type f -name '*.bak' -delete
## $ git checkout .
##
## to return FiPy back to Python 2.
