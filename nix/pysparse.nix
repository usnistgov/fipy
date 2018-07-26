{ nixpkgs, pypkgs }:
pypkgs.buildPythonPackage rec {
  pname = "pysparse";
  version = "1.1.1-dev";
  src = nixpkgs.fetchurl {
    url = "https://files.pythonhosted.org/packages/f3/a5/f1b6d0fac5bf45729e8099bd1eef039b0c7c37f8324f3be7d5e9253bb504/pysparse-1.1.1-dev.tar.gz";
    sha256 = "c291033eb8856e21b182e241f2416ee5e0565e97f2660be59ac1e3e4227e9fdb";
  };
  doCheck = false;
  buildInputs = [
    pypkgs.numpy
  ];
  hardeningDisable = [ "all" ];
}
