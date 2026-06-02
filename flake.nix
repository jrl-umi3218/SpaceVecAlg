{
  description = "SpaceVecAlg";

  inputs = {
    mc-rtc-nix.url = "github:mc-rtc/nixpkgs";
  };

  outputs =
    inputs:
    inputs.mc-rtc-nix.lib.mkFlakoboros inputs (
      { lib, ... }:
      {
        extraPackages = [ "ninja" ];
        pyOverrideAttrs.spacevecalg = {
          src = lib.cleanSource ./.;
        };
      }
    );
}
