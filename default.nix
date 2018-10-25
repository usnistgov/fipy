{ nixpkgs ? import ./nix/nixpkgs_version.nix }:
import ./nix/py2.nix { inherit nixpkgs; }
