{
  description = "SpaceVecAlg";

  inputs = {
    mc-rtc-nix.url = "github:mc-rtc/nixpkgs";
    # mc-rtc-nix.url = "path:/home/arnaud/devel/mc-rtc-nix/nixpkgs";
    # mc-rtc-nix.url = "github:arntanguy/nixpkgs-1?ref=topic/flakoboros";
  };

  outputs =
    inputs:
    inputs.mc-rtc-nix.lib.mkFlakoboros inputs (
      { lib, ... }:
      {
        extraPackages = [ "ninja" ];
        overrideAttrs.spacevecalg =
          { drv-prev, ... }:
          {
            src = lib.cleanSource ./.;
            cmakeFlags = drv-prev.cmakeFlags ++ [
              "-DPYTHON_BINDINGS=OFF"
            ];
          };
      }
    );
}
