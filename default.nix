{ stdenv
, pkgs
, cmake
, lib
, eigen
, boost
, jrl-cmakemodules
, buildPythonPackage ? null
, python313Packages ? null
, sphinx ? null
, prek ? null
, spacevecalg ? null
, isPython ? false
}:

let
  spacevecalgLib = stdenv.mkDerivation {
    pname = "spacevecalg";
    version = "0.0.0";
    src = lib.fileset.toSource {
      root = ./.;
      fileset = lib.fileset.unions [
        ./binding
        ./CMakeLists.txt
        ./src
        ./tests
        ./cmake
      ];
    };
    doCheck = true;
    nativeBuildInputs = [
      cmake
      jrl-cmakemodules
      eigen
      boost
    ];
    propagatedBuildInputs = [
      jrl-cmakemodules
    ];
    cmakeFlags = [
      "-DNANOBIND_BINDINGS=OFF"
      "-DBUILD_TESTING=ON"
    ];
  };
in
if isPython then
  buildPythonPackage rec {
    pname = "spacevecalg";
    version = "0.0.0";
    format = "other";
    src = lib.fileset.toSource {
      root = ./.;
      fileset = lib.fileset.unions [
        ./binding
        ./CMakeLists.txt
        ./src
        ./binding
        ./tests
        ./cmake
      ];
    };
    nativeBuildInputs = [
      cmake
      jrl-cmakemodules
      eigen
      boost
      python313Packages.nanobind
      python313Packages.nanoeigenpy
      sphinx
      python313Packages.sphinx-autodoc2
      python313Packages.sphinx-book-theme
      prek
      pkgs.suitesparse
    ];
    propagatedBuildInputs = [
      jrl-cmakemodules
      python313Packages.numpy
      python313Packages.nanoeigenpy
      spacevecalg
      pkgs.suitesparse
    ];
    cmakeFlags = [
      "-DNANOBIND_BINDINGS=ON"
      "-DBUILD_TESTING=OFF"
    ];
    preConfigure = ''
      export CMAKE_PREFIX_PATH="${spacevecalg}:${eigen}:${python313Packages.nanoeigenpy}:${boost}:${jrl-cmakemodules}:$CMAKE_PREFIX_PATH"
    '';
    shellHook = ''
      export PYTHONPATH="$out/${pkgs.python313.sitePackages}:$PYTHONPATH"
    '';
  }
else
  spacevecalgLib
