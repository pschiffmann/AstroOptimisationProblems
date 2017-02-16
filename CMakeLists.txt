cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
project(AstroOptimisationProblems CXX)

add_executable(gtoc1
  src/Astro_Functions.cpp
  src/Lambert.cpp
  src/mga_dsm.cpp
  src/mga.cpp
  src/misc4Tandem.cpp
  src/Pl_Eph_An.cpp
  src/PowSwingByInv.cpp
  src/propagateKEP.cpp
  src/time2distance.cpp
  src/trajobjfuns.cpp
  src/ZeroFinder.cpp

  src/main.cpp
)
set_property(TARGET gtoc1 PROPERTY CXX_STANDARD 14)
set_property(TARGET gtoc1 PROPERTY CXX_STANDARD_REQUIRED ON)
