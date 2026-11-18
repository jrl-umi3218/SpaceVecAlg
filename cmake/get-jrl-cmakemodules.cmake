# Get jrl-cmakemodules package

# Upstream (https://github.com/jrl-umi3218/jrl-cmakemodules), the new v2 version
# is located in a subfolder, We need to set this variable to bypass the v1 and
# load the v2.
set(JRL_CMAKEMODULES_USE_V2
    ON
    CACHE
      BOOL
      "Use jrl-cmakemodules v2 on https://github.com/jrl-umi3218/jrl-cmakemodules"
)

# Option 1: pass -DJRL_CMAKEMODULES_SOURCE_DIR=... to cmake command line
if(JRL_CMAKEMODULES_SOURCE_DIR)
  message(
    STATUS
      "JRL_CMAKEMODULES_SOURCE_DIR variable set, adding jrl-cmakemodules from source directory: ${JRL_CMAKEMODULES_SOURCE_DIR}"
  )
  add_subdirectory(${JRL_CMAKEMODULES_SOURCE_DIR} jrl-cmakemodules)
  return()
endif()

# Option 2: use JRL_CMAKEMODULES_SOURCE_DIR environment variable (pixi might
# unset it, prefer option 1)
if(ENV{JRL_CMAKEMODULES_SOURCE_DIR})
  message(
    STATUS
      "JRL_CMAKEMODULES_SOURCE_DIR environement variable set, adding jrl-cmakemodules from source directory: $ENV{JRL_CMAKEMODULES_SOURCE_DIR}"
  )
  add_subdirectory($ENV{JRL_CMAKEMODULES_SOURCE_DIR} jrl-cmakemodules)
  return()
endif()

# Option 3: Try to look for the installed package
message(STATUS "Looking for jrl-cmakemodules (version: >=2.1.0) package...")
find_package(jrl-cmakemodules 2.1.0 CONFIG QUIET)

# If we have the package, we are done here.
if(jrl-cmakemodules_FOUND)
  message(
    STATUS
      "Found jrl-cmakemodules (version: ${jrl-cmakemodules_VERSION}) package.")
  return()
endif()

# Option 4: Fallback to FetchContent
message(STATUS "Fetching jrl-cmakemodules using FetchContent...")
include(FetchContent)
FetchContent_Declare(
  jrl-cmakemodules
  GIT_REPOSITORY https://github.com/jrl-umi3218/jrl-cmakemodules
  GIT_TAG master)
FetchContent_MakeAvailable(jrl-cmakemodules)
