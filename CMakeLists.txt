cmake_minimum_required (VERSION 2.8)

project(IPhreeqcMMS CXX)

SET(IPhreeqcMMS_SOURCES
src/cdecl.cpp
src/fortran.cpp
src/fortran.h
src/IPhreeqcMMS.cpp
src/IPhreeqcMMS.h
src/IPhreeqcMMS.hpp
src/IPhreeqcMMSLib.cpp
src/IPhreeqcMMSLib.h
src/stdcall.cpp
)

# setup include directories
include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/IPhreeqc/src")
include_directories("${PROJECT_SOURCE_DIR}/IPhreeqc/src/phreeqcpp")
include_directories("${PROJECT_SOURCE_DIR}/IPhreeqc/src/phreeqcpp/PhreeqcKeywords")

# library
add_library(IPhreeqcMMS ${LIB_TYPE} ${IPhreeqcMMS_SOURCES})
set_target_properties(IPhreeqcMMS PROPERTIES DEBUG_POSTFIX "d")

# windows dll requires
if (MSVC AND BUILD_SHARED_LIBS)
  target_link_libraries(IPhreeqcMMS IPhreeqc)
endif()

# subdirs
add_subdirectory(IPhreeqc)
