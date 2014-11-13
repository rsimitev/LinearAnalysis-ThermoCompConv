# CMake project file for the Glasgow linear onset code.

cmake_minimum_required (VERSION 2.6)
project(GLO)
enable_language (Fortran)
set (GLO_VERSION_MAJOR 0)
set (GLO_VERSION_MINOR 6)
set (GLO_VERSION_PATCH 0)
set (GLO_VERSION
   "${GLO_VERSION_MAJOR}.${GLO_VERSION_MINOR}.${GLO_VERSION_PATCH}" )

SET(CMAKE_MODULE_PATH ${GLO_SOURCE_DIR}/cmake)
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-cpp -ffree-line-length-none -fno-f2c -O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-cpp -ffree-line-length-none -fno-f2c -O0 -g -Wall -Wextra -fcheck=all -finit-real=inf")
  set (CMAKE_Fortran_FLAGS_PROFILE "-cpp -ffree-line-length-none -fno-f2c -O0 -g -pg -Wall -Wextra")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2 -fpp -free -heap-arrays")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fpp -free -warn all -check bounds -traceback  -fp-stack-check -heap-arrays")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")

set (GLO_VERSION_STRING "${GLO_VERSION}(${FEATURES}) - ${CMAKE_BUILD_TYPE}")

include_directories(${GLO_BINARY_DIR})

set(BLA_STATIC ON)
find_package(LAPACK REQUIRED)

add_subdirectory(src)
add_subdirectory(scripts)

# add a target to generate API documentation with Doxygen
find_package(Doxygen)
if(DOXYGEN_FOUND)
   configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
      ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
   add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
      )
endif(DOXYGEN_FOUND)

INCLUDE(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/BUILD/;/bin/;test*;/.git*;/.kate*;run-test*")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${GLO_SOURCE_DIR}/README") 
set(CPACK_PACKAGE_VENDOR "The University of Glasgow")
set(CPACK_PACKAGE_VERSION_MAJOR ${GLO_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${GLO_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${GLO_VERSION_PATCH})

INCLUDE(CPack)

