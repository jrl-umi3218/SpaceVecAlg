{
  lib,

  stdenv,

  # nativeBuildInputs
  cmake,
  eigen,
  boost,
  python3Packages,

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
    ];
  };

  doCecks = true; # run unit tests
  nativeBuildInputs = [
    cmake
    eigen
    boost
    python3Packages.python
  ];

  propagatedBuildInputs = [
    jrl-cmakemodules
    python3Packages.numpy
  ];

  # FIXME: upgrade python bindings (nanobind?)
  cmakeFlags = [
    "-DPYTHON_BINDING=OFF"
    "-DBUILD_TESTING=ON"
  ];
}
