{
  description = "Implementation of spatial vector algebra with the Eigen3 linear algebra library.";

  # inputs.mc-rtc-nix.url = "github:mc-rtc/nixpkgs";
  inputs.mc-rtc-nix.url = "github:mc-rtc/nixpkgs/pull/33/head";
  # inputs.mc-rtc-nix.url = "path:/home/arnaud/devel/mc-rtc-nix/nixpkgs";

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
            pname = "spacevecalg-nanobind";
            src = lib.cleanSource ./.;
            outputs = [
              "out"
              "doc"
            ];
            cmakeFlags = [
              (lib.cmakeBool "PYTHON_BINDINGS" false)
              (lib.cmakeBool "BUILD_DOCUMENTATION" true)
              (lib.cmakeBool "INSTALL_DOCUMENTATION" true)
              (lib.cmakeBool "NANOBIND_BINDINGS" true)
              (lib.cmakeBool "NANOBIND_DOCUMENTATION" true)
              (lib.cmakeBool "BUILD_TESTING" false)
            ];
            nativeBuildInputs = with pkgs-final; [
              cmake
              jrl-cmakemodulesv2
              doxygen
              graphviz
              python3Packages.python
              # nanobind documentation
              sphinx
              sphinx-cmake
              python3Packages.sphinx-autodoc2
              python3Packages.sphinx-book-theme
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
