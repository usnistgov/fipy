let
    pkgs = import (builtins.fetchGit {
      url = "https://github.com/NixOS/nixpkgs.git";
      rev = "92e1376cc3e37ee72f1417df05032556d39853c1";
      ref = "master";
    }) { };
    pythonPackages = pkgs.python3Packages;
    not_darwin_inputs = pkgs.lib.optionals (! pkgs.stdenv.isDarwin ) [ pythonPackages.jupyter ];
    not_darwin_pre_shell_hook = if (! pkgs.stdenv.isDarwin) then ''
      jupyter nbextension install --py widgetsnbextension --user
      jupyter nbextension enable widgetsnbextension --user --py
    '' else "";
in
  (pythonPackages.fipy.overridePythonAttrs (old: rec {
    src = builtins.filterSource (path: type: type != "directory" || baseNameOf path != ".git") ./.;
    nativeBuildInputs = with pythonPackages; [
      pip
      pkgs.imagemagick
      pkgs.git
    ] ++ old.propagatedBuildInputs ++ not_darwin_inputs;
    postShellHook = not_darwin_pre_shell_hook + ''
      SOURCE_DATE_EPOCH=$(date +%s)
      export PYTHONUSERBASE=$PWD/.local
      export USER_SITE=`python -c "import site; print(site.USER_SITE)"`
      export PYTHONPATH=$PYTHONPATH:$USER_SITE
      export PATH=$PATH:$PYTHONUSERBASE/bin

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
