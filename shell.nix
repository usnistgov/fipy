let
    pkgs = import (builtins.fetchGit {
      url = "https://github.com/NixOS/nixpkgs.git";
      rev = "c2ae05d5973cc4f8842755f2807ac10e31bb2aa8";
      ref = "master";
    }) { };
    pythonPackages = pkgs.python3Packages;
    not_darwin_inputs = pkgs.lib.optionals (! pkgs.stdenv.isDarwin ) [ pythonPackages.jupyter ];
    not_darwin_pre_shell_hook = if (! pkgs.stdenv.isDarwin) then ''
      jupyter nbextension install --py widgetsnbextension --user
      jupyter nbextension enable widgetsnbextension --user --py
    '' else "";
    filter_pyamg = builtins.filter (x: ! ((pkgs.lib.hasInfix "pyamg" x.name) && pkgs.stdenv.isDarwin));
in
  (pythonPackages.fipy.overridePythonAttrs (old: rec {
    src = builtins.filterSource (path: type: type != "directory" || baseNameOf path != ".git") ./.;
    nativeBuildInputs = with pythonPackages; [
      pip
      pkgs.imagemagick
      pkgs.git
    ] ++ propagatedBuildInputs ++ not_darwin_inputs;

    propagatedBuildInputs = (filter_pyamg old.propagatedBuildInputs);

    postShellHook = not_darwin_pre_shell_hook + ''
      SOURCE_DATE_EPOCH=$(date +%s)
      export PYTHONUSERBASE=$PWD/.local
      export USER_SITE=`python -c "import site; print(site.USER_SITE)"`
      export PYTHONPATH=$PYTHONPATH:$USER_SITE
      export PATH=$PATH:$PYTHONUSERBASE/bin

      export OMPI_MCA_plm_rsh_agent=/usr/bin/ssh

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
