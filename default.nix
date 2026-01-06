{
  lib,

  stdenv,

  # nativeBuildInputs
  cmake,
  eigen,
  boost,
  git,
  python313Packages,

  # propagatedBuildInputs
  jrl-cmakemodules,
}:
stdenv.mkDerivation {
  pname = "spacevecalg";
  version = "0.0.0";

  src = lib.fileset.toSource {
    root = ./.;
    fileset = lib.fileset.unions [
        ./binding
        ./CMakeLists.txt
        ./src
        ./tests
        ./cmake-v2
    ];
  };

  doCecks = true; # run unit tests
  nativeBuildInputs = [
    cmake
    eigen
    boost
    git
    python313Packages.python
    python313Packages.nanobind
  ];

  propagatedBuildInputs = [
    jrl-cmakemodules
    python313Packages.numpy
  ];

  # FIXME: upgrade python bindings (nanobind?)
  cmakeFlags = [
    "-DPYTHON_BINDING=OFF"
    "-DNANOBIND_BINDINGS=ON"
    "-DBUILD_TESTING=ON"
  ];
}
