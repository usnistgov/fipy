{ nixpkgs ? import ./nixpkgs_version.nix }:
let
  pypkgs = nixpkgs.python27Packages;
  pysparse = import ./pysparse.nix { inherit nixpkgs pypkgs; };
  preshellhook = _: '' '';
in
  import ./build.nix { inherit nixpkgs pypkgs pysparse preshellhook; }
