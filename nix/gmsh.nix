{ nixpkgs }:
let version = "3.0.6"; in
nixpkgs.stdenv.mkDerivation {
  name = "gmsh-${version}";

  src = nixpkgs.fetchurl {
    url = "http://gmsh.info/src/gmsh-${version}-source.tgz";
    sha256 = "0ywqhr0zmdhn8dvi6l8z1vkfycyv67fdrz6b95mb39np832bq04p";
  };

  # The original CMakeLists tries to use some version of the Lapack lib
  # that is supposed to work without Fortran but didn't for me.
  patches = [ ./CMakeLists.txt.patch ];

  buildInputs = [
    nixpkgs.pkgs.cmake
    nixpkgs.pkgs.blas
    nixpkgs.pkgs.liblapack
    nixpkgs.pkgs.gfortran
    nixpkgs.pkgs.gmm
    nixpkgs.pkgs.fltk
    nixpkgs.pkgs.libjpeg
    nixpkgs.pkgs.zlib
    nixpkgs.pkgs.libGLU_combined
    nixpkgs.pkgs.libGLU
    nixpkgs.pkgs.xorg.libXrender
    nixpkgs.pkgs.xorg.libXcursor
    nixpkgs.pkgs.xorg.libXfixes
    nixpkgs.pkgs.xorg.libXext
    nixpkgs.pkgs.xorg.libXft
    nixpkgs.pkgs.xorg.libXinerama
    nixpkgs.pkgs.xorg.libX11
    nixpkgs.pkgs.xorg.libSM
    nixpkgs.pkgs.xorg.libICE
    nixpkgs.gfortran.cc.lib
  ];

  enableParallelBuilding = true;

  meta = {
    description = "A three-dimensional finite element mesh generator";
    homepage = http://gmsh.info/;
    platforms = nixpkgs.stdenv.lib.platforms.all;
    license = nixpkgs.stdenv.lib.licenses.gpl2Plus;
  };
}
