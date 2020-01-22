cmake_minimum_required (VERSION 2.8.11)

project(webmod Fortran C CXX)

SET(webmod_SOURCES
../IPhreeqcMMS/include/mms_phreeqc.f90.inc
../IPhreeqcMMS/IPhreeqc/src/IPhreeqc_interface.F90
../IPhreeqcMMS/src/phr_mix.F90
../IPhreeqcMMS/src/phr_multicopy.f90
../IPhreeqcMMS/src/phr_precip.f90
src/basin_topg.F
src/call_modules_webmod.F
src/ccsolrad_web.F
src/defines.h
src/elapsed_time_function.F90
src/fmodules.inc
src/fmodules.f90.inc
src/intcp_prms.F
src/io.F
src/irrig_web.F
src/my_newunit.F
src/nwsmelt_topg.F
src/obs_chem.F
src/obs_webmod.F
src/phreeq_mms.F90
src/potet_hamon_prms.F
src/precip_web.F
src/route_clark.F
src/soltab_prms.F
src/temp_1sta_prms.F
src/top2clark.F
src/topmod_chem.F
src/web_sum.F
src/webmod_res.F
)

# define CMAKE_FC
add_definitions(-DCMAKE_FC)

# place [*.F[90]] files into "Source Files"
source_group("Header Files" FILES "../IPhreeqcMMS/include/mms_phreeqc.f90.inc")
source_group("Source Files" FILES "../IPhreeqcMMS/IPhreeqc/src/IPhreeqc_interface.F90")
source_group("Source Files" FILES "../IPhreeqcMMS/src/phr_mix.F90")
source_group("Source Files" FILES "src/basin_topg.F")
source_group("Source Files" FILES "src/call_modules_webmod.F")
source_group("Source Files" FILES "src/ccsolrad_web.F")
source_group("Header Files" FILES "src/defines.h")
source_group("Source Files" FILES "src/elapsed_time_function.F90")
source_group("Header Files" FILES "src/fmodules.inc")
source_group("Header Files" FILES "src/fmodules.f90.inc")
source_group("Source Files" FILES "src/intcp_prms.F")
source_group("Source Files" FILES "src/io.F")
source_group("Source Files" FILES "src/irrig_web.F")
source_group("Source Files" FILES "src/my_newunit.F")
source_group("Source Files" FILES "src/nwsmelt_topg.F")
source_group("Source Files" FILES "src/obs_chem.F")
source_group("Source Files" FILES "src/obs_webmod.F")
source_group("Source Files" FILES "src/phreeq_mms.F90")
source_group("Source Files" FILES "src/potet_hamon_prms.F")
source_group("Source Files" FILES "src/precip_web.F")
source_group("Source Files" FILES "src/route_clark.F")
source_group("Source Files" FILES "src/soltab_prms.F")
source_group("Source Files" FILES "src/temp_1sta_prms.F")
source_group("Source Files" FILES "src/top2clark.F")
source_group("Source Files" FILES "src/topmod_chem.F")
source_group("Source Files" FILES "src/web_sum.F")
source_group("Source Files" FILES "src/webmod_res.F")

# setup include directories
include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/../IPhreeqcMMS/include")
include_directories("${PROJECT_SOURCE_DIR}/../IPhreeqcMMS/IPhreeqc/src")

# make static
if (MSVC)
  set(CompilerFlags
    CMAKE_CXX_FLAGS
    CMAKE_CXX_FLAGS_DEBUG
    CMAKE_CXX_FLAGS_RELEASE
    CMAKE_CXX_FLAGS_MINSIZEREL
    CMAKE_CXX_FLAGS_RELWITHDEBINFO
    CMAKE_C_FLAGS
    CMAKE_C_FLAGS_DEBUG
    CMAKE_C_FLAGS_RELEASE
    CMAKE_C_FLAGS_MINSIZEREL
    CMAKE_C_FLAGS_RELWITHDEBINFO
  )

  # Set C/C++ Runtime Library to Multi-threaded (/MT) 
  foreach(CompilerFlag ${CompilerFlags})
    if(${CompilerFlag} MATCHES "/MD")
      string(REGEX REPLACE "/MD" "/MT" ${CompilerFlag} "${${CompilerFlag}}")
    endif()
  endforeach()
  
  if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} /MP")
    
    set(FortranCompilerFlags
      CMAKE_Fortran_FLAGS
    )
  
    # Set Fortran Runtime Library to Multithreaded
    foreach(FortranCompilerFlag ${FortranCompilerFlags})
      if(${FortranCompilerFlag} MATCHES "/libs:dll")
        string(REGEX REPLACE "/libs:dll" "/libs:static" ${FortranCompilerFlag} "${${FortranCompilerFlag}}")
      endif()  
    endforeach()
  endif()
endif()

# executable
add_executable(webmod ${webmod_SOURCES})
set_target_properties(webmod PROPERTIES DEBUG_POSTFIX "d")
###set_target_properties(webmod PROPERTIES LINKER_LANGUAGE "CXX")

# library dependencies
set(EXTRA_LIBS ${EXTRA_LIBS} IPhreeqcMMS)
set(EXTRA_LIBS ${EXTRA_LIBS} IPhreeqc)
set(EXTRA_LIBS ${EXTRA_LIBS} mmf_c)

# link 
target_link_libraries(webmod ${EXTRA_LIBS})

# install
install(TARGETS webmod DESTINATION ${CMAKE_INSTALL_BINDIR})

