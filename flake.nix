{
  description = "Implementation of spatial vector algebra with the Eigen3 linear algebra library.";

  inputs.mc-rtc-nix.url = "github:mc-rtc/nixpkgs";

  outputs =
    inputs:
    inputs.mc-rtc-nix.lib.mkFlakoboros inputs (
      { lib, ... }:
      {
        extraDevPyPackages = [ "spacevecalg" ];

        overrides.spacevecalg =
          { pkgs-final, ... }:
          {
            jrl-cmakemodules = pkgs-final.jrl-cmakemodulesv2;
          };
        overrideAttrs.spacevecalg =
          { pkgs-final, ... }:
          {
            pname = "spacevecalg";
            src = lib.cleanSource ./.;
            cmakeFlags = [
              (lib.cmakeBool "PYTHON_BINDINGS" false)
              (lib.cmakeBool "NANOBIND_BINDINGS" true)
              (lib.cmakeBool "BUILD_TESTING" false)
            ];
            nativeBuildInputs = with pkgs-final; [
              cmake
              jrl-cmakemodulesv2
              doxygen
              python3Packages.python
            ];
            propagatedBuildInputs = with pkgs-final; [
              eigen
              boost
              python3Packages.nanoeigenpy
              python3Packages.nanobind
            ];
          };

        pyPackages = {
          spacevecalg =
            {
              pkgs,
              toPythonModule,
            }:
            (toPythonModule (pkgs.spacevecalg.override { }));
        };
      }
    );
}
