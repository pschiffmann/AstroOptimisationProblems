cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(AstroOptimisationProblems CXX)

find_path(MANTELLA_INCLUDE_DIR mantella0)

add_executable(gtoc1
  src/Astro_Functions.cpp
  src/ZeroFinder.cpp

  src/mantella_bindings/gtoc1.cpp
  src/mantella_bindings/astro_helpers.cpp
  src/mantella_bindings/vector3d_helpers.cpp
  src/mantella_bindings/main.cpp
)
set_property(TARGET gtoc1 PROPERTY CXX_STANDARD 14)
set_property(TARGET gtoc1 PROPERTY CXX_STANDARD_REQUIRED ON)
target_include_directories(gtoc1 SYSTEM PUBLIC ${MANTELLA_INCLUDE_DIR})
