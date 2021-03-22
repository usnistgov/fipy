#
# $ nix-shell --pure --argstr tag 20.09
#
{
  tag ? "20.09"
}:
let
  pkgs = import (builtins.fetchTarball "https://github.com/NixOS/nixpkgs/archive/${tag}.tar.gz") {};
  pypkgs = pkgs.python3Packages;
  not_darwin_inputs = pkgs.lib.optionals (! pkgs.stdenv.isDarwin ) [ pypkgs.jupyter ];
  not_darwin_pre_shell_hook = if (! pkgs.stdenv.isDarwin) then ''
    jupyter nbextension install --py widgetsnbextension --user
    jupyter nbextension enable widgetsnbextension --user --py
  '' else "";
  filter_pyamg = builtins.filter (x: ! ((pkgs.lib.hasInfix "pyamg" x.name) && pkgs.stdenv.isDarwin));
in
  (pypkgs.fipy.overridePythonAttrs (old: rec {
    src = builtins.filterSource (path: type: type != "directory" || baseNameOf path != ".git") ./.;
    nativeBuildInputs = with pypkgs; [
      pip
      pkgs.imagemagick
      pkgs.git
      pkgs.openssh
    ] ++ propagatedBuildInputs ++ not_darwin_inputs;

    propagatedBuildInputs = (filter_pyamg old.propagatedBuildInputs);

    postShellHook = not_darwin_pre_shell_hook + ''
      SOURCE_DATE_EPOCH=$(date +%s)
      export PYTHONUSERBASE=$PWD/.local
      export USER_SITE=`python -c "import site; print(site.USER_SITE)"`
      export PYTHONPATH=$PYTHONPATH:$USER_SITE
      export PATH=$PATH:$PYTHONUSERBASE/bin

      export OMPI_MCA_plm_rsh_agent=${pkgs.openssh}/bin/ssh

      ## To build the docs
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
