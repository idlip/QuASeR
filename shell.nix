with import <nixpkgs> {};
pkgs.mkShell {

  nativeBuildInputs = [ pkgs.bashInteractive ];

  # EnvVars = The thung
  # NIX_LD_LIBRARY_PATH = lib.makeLibraryPath [
  # pkgs
  # ];
  # NIX_LD = lib.fileContents "${stdenv.cc}/nix-support/dynamic-linker";

  buildInputs = with pkgs; [
    uv git
  ];

  shellHook = ''
    echo "You've entered Quaser development Dungeon"
  '';
}
