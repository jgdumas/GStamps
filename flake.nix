{
  description = "Flake GStamps project (https://github.com/jgdumas/GStamps/).";

  # To easily generate a derivation per architecture
  inputs.flake-utils.url = "github:numtide/flake-utils";
  
  outputs = { self, nixpkgs, flake-utils }: {
    # We define first an overlay, i.e. a definition of new packages as recommended in
    # https://discourse.nixos.org/t/how-to-consume-a-eachdefaultsystem-flake-overlay/19420/9
    overlays.default = final: prev: {
      gstamps = final.callPackage ({stdenv, pkgs}:
        stdenv.mkDerivation {
          src = ./.;
          pname = "gstamps";
          version = "unstable";
          # Warning: this introduces some impurity as it depends on the system it is running on (optimization -march=native).
          # For this reason you may want to disable this line, and/or change the Makefile to compile for a precise architecture.
          NIX_ENFORCE_NO_NATIVE = false;
          buildInputs = with pkgs; [
            boost
            pkg-config
            givaro
          ];
          # The buildPhase is automatically running "make all", so we just need to write the installPhase:
          installPhase = ''
            mkdir -p $out/bin
            cp -r bin/* $out/bin
          '';
        }
      ) {};
    };
  } // # We now add the package defined in the above overlay for all architectures
  (flake-utils.lib.eachDefaultSystem (system:
    let
      pkgs = import nixpkgs {
        system = system;
        overlays = [ self.overlays.default ];
      };
      lib = nixpkgs.lib;
    in
      {
        # Create a new package
        packages = {
          default = pkgs.gstamps;
          gstamps = pkgs.gstamps;
        };
      }
  ))
  ;
}
