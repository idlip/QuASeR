with import <nixpkgs> {};
pkgs.mkShell {

  nativeBuildInputs = [ pkgs.bashInteractive ];

  buildInputs = with pkgs; [
    uv git
  ];

  shellHook = ''
    echo "You've entered Quaser development Dungeon."
  '';
}
