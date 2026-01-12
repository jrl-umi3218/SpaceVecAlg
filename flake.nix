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
      systems = import inputs.systems;
      imports = [ inputs.gepetto.flakeModule ];
      perSystem =
        {
          pkgs,
          self',
          inputs',
          ...
        }:
        {
          packages = {
            default = self'.packages.spacevecalg;
            spacevecalg = pkgs.callPackage ./. {
                jrl-cmakemodules = inputs'.jrl-cmakemodules.packages.default;
            };
          };
        };
    };
}
