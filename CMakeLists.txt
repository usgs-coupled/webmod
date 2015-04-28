cmake_minimum_required (VERSION 2.8.11)

project(web_win C CXX Fortran)

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
# enable_testing()

# subdirs
add_subdirectory(IPhreeqcMMS)
add_subdirectory(mmf_c)
add_subdirectory(webmod)
