{ pypkgs }:
pypkgs.buildPythonPackage rec {
  version = "0.0.9";
  pname = "scikit-fmm";
  src = pypkgs.fetchPypi {
    inherit pname version;
    sha256 = "192f4rzh32bz3gs5cqybxrg0jqjyp751kxfdx27paf81v7cwg2ky";
  };
  doCheck=false;
  buildInputs = [
    pypkgs.numpy
  ];
}
