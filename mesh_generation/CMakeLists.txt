cmake_minimum_required(VERSION 3.0.0 FATAL_ERROR)

project(ccs LANGUAGES Fortran)

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -r8 -parallel")
  if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -gdwarf -Og -traceback -check all")
  endif ()
elseif (CMAKE_Fortran_COMPILER_ID MATCHES GNU)
  if (APPLE)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -floop-parallelize-all -fdefault-real-8 -ffree-line-length-none")
  else ()
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O3 -floop-parallelize-all -fdefault-real-8 -ftree-parallelize-loops=8 -ffree-line-length-none -lgomp -fopenmp")
  endif ()
  if (CMAKE_BUILD_TYPE MATCHES DEBUG)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -gdwarf -Og -fbacktrace -fcheck=all")
  endif ()
endif ()

if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()
if (DEFINED ENV{NETCDF_ROOT})
  include_directories("$ENV{NETCDF_ROOT}/include")
  link_directories("$ENV{NETCDF_ROOT}/lib")
else ()
  find_package(PkgConfig)
  if (PKG_CONFIG_FOUND)
    pkg_search_module(NETCDF REQUIRED netcdf)
    if (NETCDF_FOUND)
      include_directories(${NETCDF_INCLUDEDIR})
      link_directories(${NETCDF_LIBDIR})
    else ()
      message(FATAL_ERROR "Pkg-config could not find netcdf library!")
    endif ()
  else ()
    message(FATAL_ERROR "Unable to find pkg-config library!")
  endif ()
endif ()

set(source_files
  src/ccs.f90
  src/constants_mod.f90
  src/mesh_mod.f90
  src/math_mod.f90
  src/ghost_mod.f90
  src/output_mod.f90
  src/parameters_mod.f90
  src/cuco.f90
  src/pcoco.f90
  src/pfft.f90
  src/pfft1.f90
  src/pietc.f90
  src/pietc_s.f90
  src/pkind.f90
  src/pmat.f90
  src/pmat4.f90
  src/pmat5.f90
  src/jimcc_m.f90
  src/precis_m.f90
  src/parm_m.f90
)

add_library(ccs ${source_files})
target_link_libraries(ccs netcdf netcdff)

add_executable(ccs.exe src/ccs.F90)
target_link_libraries(ccs.exe ccs)

