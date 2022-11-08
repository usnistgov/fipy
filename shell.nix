#
# $ nix-shell --pure --argstr tag 20.09
#
{
  tag ? "22.05"
}:
let
  pkgs = import (builtins.fetchTarball "https://github.com/NixOS/nixpkgs/archive/${tag}.tar.gz") {};
  pypkgs = pkgs.python3Packages;

  nixes_src = builtins.fetchTarball "https://github.com/wd15/nixes/archive/9a757526887dfd56c6665290b902f93c422fd6b1.zip";
  jupyter_extra = pypkgs.callPackage "${nixes_src}/jupyter/default.nix" {
    jupyterlab=(if pkgs.stdenv.isDarwin then pypkgs.jupyter else pypkgs.jupyterlab);
  };

in
  (pypkgs.fipy.overridePythonAttrs (old: rec {

    src = pkgs.lib.cleanSource ./.;

    nativeBuildInputs = with pypkgs; [
      pip
      pkgs.imagemagick
      pkgs.git
      pkgs.openssh
      nbval
    ] ++ propagatedBuildInputs ++ [ jupyter_extra ipywidgets ];

    propagatedBuildInputs = old.propagatedBuildInputs;

    postShellHook = ''
      SOURCE_DATE_EPOCH=$(date +%s)
      export PYTHONUSERBASE=$PWD/.local
      export USER_SITE=`python -c "import site; print(site.USER_SITE)"`
      export PYTHONPATH=$PYTHONPATH:$USER_SITE:$(pwd)
      export PATH=$PATH:$PYTHONUSERBASE/bin

      export OMPI_MCA_plm_rsh_agent=${pkgs.openssh}/bin/ssh

      ## To build the docs (don't use --pure)
      # Need latex to build properly
      #
      # texlive-latex-base, texlive-fonts-recommended, texlive-fonts-extra
      # texlive-latex-extra, texlive-science, texlive-extra-utils
      #
      # pip install --user sphinx
      # pip install --user "sphinxcontrib-bibtex<=0.4.2"
      # pip install --user git+https://github.com/thewtex/sphinx-contrib.git#subdirectory=traclinks
      # required for embedded plots in documentation
      # pip install --user pandas
      # python setup.py build_docs --html

      ## To build PyAMG add nixpkgs.gcc to nativeBuildInputs and then
      # pip install --user pyamg
    '';
  }))
