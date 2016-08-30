cmake_minimum_required (VERSION 2.8.11)

project(webmod C CXX Fortran)

# set version numbers
set(WEBMOD_VERSION_MAJOR "1")
set(WEBMOD_VERSION_MINOR "0")
set(WEBMOD_VERSION_PATCH "0")
set(WEBMOD_REVISION "9000")

# set default windows installation path
if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT AND WIN32)
  GetDefaultWindowsPrefixBase(CMAKE_GENERIC_PROGRAM_FILES)  
  SET(CMAKE_INSTALL_PREFIX "${CMAKE_GENERIC_PROGRAM_FILES}/USGS/${PROJECT_NAME}-${WEBMOD_VERSION_MAJOR}.${WEBMOD_VERSION_MINOR}.${WEBMOD_VERSION_PATCH}}" CACHE PATH "Install path prefix, prepeneded onto install directories." FORCE)   
endif()

# overide docdir on windows
if(WIN32 AND NOT CMAKE_INSTALL_DOCDIR)
  set(CMAKE_INSTALL_DOCDIR "" CACHE PATH "documentation root (doc)")
  set(CMAKE_INSTALL_DOCDIR "doc")
endif()

# retrieve default install locations
include(GNUInstallDirs)

# Disable shared libraries on windows
# see https://support2.microsoft.com/default.aspx?scid=KB;EN-US;168958
# and https://msdn.microsoft.com/query/dev11.query?appId=Dev11IDEF1&l=EN-US&k=k(C4251);&rd=true
if (WIN32)
  set(LIB_TYPE STATIC)
else()
  if(NOT DEFINED BUILD_SHARED_LIBS)
    option(BUILD_SHARED_LIBS "Build Shared Libraries" OFF)
    set(LIB_TYPE STATIC)
  endif()
endif()

if (BUILD_SHARED_LIBS)
  set(LIB_TYPE SHARED)
endif()

# Disable deprecation warnings for standard C and STL functions in VS2005
# and later
if(MSVC_VERSION EQUAL 1400 OR MSVC_VERSION GREATER 1400)
  add_definitions(-D_CRT_SECURE_NO_DEPRECATE)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
endif(MSVC_VERSION EQUAL 1400 OR MSVC_VERSION GREATER 1400) 

# tests
enable_testing()

# subdirs
add_subdirectory(IPhreeqcMMS)
add_subdirectory(mmf_c)
add_subdirectory(webmod.d)
add_subdirectory(Tests)
add_subdirectory(doc)
add_subdirectory(lib)
