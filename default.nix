{
  lib,

  stdenv,

  # nativeBuildInputs
  cmake,
  eigen,
  boost,
  python313Packages,
  sphinx,

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
        ./cmake
    ];
  };

  doCecks = true; # run unit tests
  nativeBuildInputs = [
    # build dependencies
    cmake
    jrl-cmakemodules
    eigen
    boost
    python313Packages.python
    python313Packages.nanobind
    python313Packages.nanoeigenpy
    # testing
    python313Packages.pytest
    python313Packages.scipy
    # documentation
    sphinx
    python313Packages.sphinx-autodoc2
    python313Packages.sphinx-book-theme
  ];

  propagatedBuildInputs = [
    jrl-cmakemodules
    python313Packages.numpy
  ];

  # FIXME: upgrade python bindings (nanobind?)
  cmakeFlags = [
    "-DNANOBIND_BINDINGS=ON"
    "-DBUILD_TESTING=ON"
  ];

  shellHook = ''
    export PYTHONPATH="$PWD/build/lib/site-packages:$PYTHONPATH"
  '';
}
