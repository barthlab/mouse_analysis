{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs";
  outputs = {
    self,
    nixpkgs,
  }: {
    devShells.x86_64-linux.default = with import nixpkgs {
      system = "x86_64-linux";
      config = {
        allowUnfree = true;
      };
    };
      mkShell {
        buildInputs = [
          pkgs.python3Full
          pkgs.python3Packages.matplotlib
          pkgs.python3Packages.pandas
        ];

        shellHook = ''
          ${pkgs.zsh}/bin/zsh
          exit'';
      };
  };
}
