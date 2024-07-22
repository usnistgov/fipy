{
  description = "Python environment for fipy";

  inputs = {
    nixpkgs.url = "github:NixOS/nixpkgs/de02e640d0e7787441326133b6de1b44b2d3865f";
    utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, utils}: (utils.lib.eachSystem [ "x86_64-linux" "x86_64-darwin" "aarch64-darwin" ] (system:
    let
      pkgs = nixpkgs.legacyPackages.${system};
      pypkgs = pkgs.python312Packages;

      env = (pypkgs.fipy.overridePythonAttrs (old: rec {

        src = pkgs.lib.cleanSource ./.;

        disabled = pypkgs.pythonOlder "3.7";

        nativeBuildInputs_ = with pypkgs; [
          pip
          pkgs.openssh
          nbval
          ipython
          ipykernel
          jupyterlab
          traitlets
          notebook
          scipy
        ];

        nativeBuildInputs = propagatedBuildInputs;

        propagatedBuildInputs = old.propagatedBuildInputs ++ nativeBuildInputs_;

        postShellHook = ''
          SOURCE_DATE_EPOCH=$(date +%s)
          export PYTHONUSERBASE=$PWD/.local
          export USER_SITE=`python -c "import site; print(site.USER_SITE)"`
          export PYTHONPATH=$PYTHONPATH:$USER_SITE:$(pwd):$PWD
          export PATH=$PATH:$PYTHONUSERBASE/bin

          export OMPI_MCA_plm_rsh_agent=${pkgs.openssh}/bin/ssh
        '';
      }));
   in
     rec {
       packages.fipy = env;
       packages.default = self.packages.${system}.fipy;
       devShells.default = pkgs.mkShell {
        packages = [ env ];
       };
       devShells.develop = env;
     }
  ));
}
