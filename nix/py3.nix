{ nixpkgs ? import ./nixpkgs_version.nix }:
let
  pypkgs = nixpkgs.python36Packages;
  pysparse = null;
  preshellhook = _: '' '';
in
  import ./build.nix { inherit nixpkgs pypkgs pysparse preshellhook; }
