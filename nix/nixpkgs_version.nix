let
  inherit (import <nixpkgs> {}) fetchFromGitHub;
  nixpkgs_download = fetchFromGitHub {
    owner = "NixOS";
    repo = "nixpkgs-channels";
    rev = "120b013e0c082d58a5712cde0a7371ae8b25a601";
    sha256 = "00gd96p7yz3rgpjjkizp397y2syfc272yvwxqixbjd1qdshbizmj";
  };
in
  import nixpkgs_download {}
