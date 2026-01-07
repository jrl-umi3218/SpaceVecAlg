# FIXME: aarch64, python bindings
{
  description = "Implementation of spatial vector algebra with the Eigen3 linear algebra library.";

  inputs = {
    gepetto.url = "github:gepetto/nix";
    flake-parts.follows = "gepetto/flake-parts";
    nixpkgs.follows = "gepetto/nixpkgs";
    nix-ros-overlay.follows = "gepetto/nix-ros-overlay";
    systems.follows = "gepetto/systems";
    treefmt-nix.follows = "gepetto/treefmt-nix";
    # XXX: use official flake when v2 PR is merged
    jrl-cmakemodules = {
      url = "github:ahoarau/jrl-cmakemodules?ref=jrl-next";
    };
  };

  outputs =
    inputs:
    inputs.flake-parts.lib.mkFlake { inherit inputs; } {
      # FIXME: cross-compilation for aarch64 does not work as-is
      # systems = import inputs.systems;
      systems = [ "x86_64-linux" ];
      imports = [ inputs.gepetto.flakeModule ];
      perSystem =
        {
          pkgs,
          self',
          ...
        }:
        {
          packages = {
            default = self'.packages.spacevecalg;
            spacevecalg = pkgs.callPackage ./. {
                jrl-cmakemodules = inputs.jrl-cmakemodules.packages.${pkgs.system}.default;
            };
          };
        };
    };
}
