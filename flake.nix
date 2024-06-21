{
  description = "Python environment for fipy";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.05";
    utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, utils}: (utils.lib.eachSystem [ "x86_64-linux" "x86_64-darwin" "aarch64-darwin" ] (system:
    let
      pkgs = nixpkgs.legacyPackages.${system};
      pypkgs = pkgs.python312Packages;

      env = (pypkgs.fipy.overridePythonAttrs (old: rec {

        src = pkgs.lib.cleanSource ./.;

        nativeBuildInputs = with pypkgs; [
          pip
          pkgs.openssh
          nbval
          ipython
          ipykernel
          jupyterlab
          traitlets
          notebook
        ] ++ propagatedBuildInputs;

        propagatedBuildInputs = old.propagatedBuildInputs;

        postShellHook = ''
          SOURCE_DATE_EPOCH=$(date +%s)
          export PYTHONUSERBASE=$PWD/.local
          export USER_SITE=`python -c "import site; print(site.USER_SITE)"`
          export PYTHONPATH=$PYTHONPATH:$USER_SITE:$(pwd)
          export PATH=$PATH:$PYTHONUSERBASE/bin

          export OMPI_MCA_plm_rsh_agent=${pkgs.openssh}/bin/ssh

        '';
      }));
   in
     rec {
       packages.fipy = env;
       packages.default = self.packages.${system}.fipy;
       devShells.default = env;
     }
  ));
}
