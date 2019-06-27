let
  inherit (import <nixpkgs> {}) fetchFromGitHub;
  nixpkgs_download = fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs-channels";
    # rev = "61deecdc34fc609d0f805b434101f3c8ae3b807a";
    # sha256 = "147xyn8brvkfgz1z3jbk13w00h6camnf6z0bz0r21g9pn1vv7sb0";
    rev = "f52505fac8c82716872a616c501ad9eff188f97f"; # version: 19.03
    sha256 = "0q2m2qhyga9yq29yz90ywgjbn9hdahs7i8wwlq7b55rdbyiwa5dy";
  };
in
  import nixpkgs_download {}
