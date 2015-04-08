cmake_minimum_required (VERSION 2.8)

project(webmod_f Fortran)

##include(FortranCInterface)
##FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_")
##FortranCInterface_VERIFY(CXX QUIET)

SET(webmod_f_SOURCES
../IPhreeqcMMS/include/mms_phreeqc.f90.inc
../IPhreeqcMMS/IPhreeqc/src/IPhreeqc_interface.F90
../IPhreeqcMMS/src/phr_mix.F90
../IPhreeqcMMS/src/phr_multicopy.f90
../IPhreeqcMMS/src/phr_precip.f90
src/basin_topg.f
src/call_modules_webmod.f
src/ccsolrad_web.f
src/elapsed_time_function.f90
src/fmodules.f90.inc
src/intcp_prms.f
src/io.f
src/irrig_web.f
src/nwsmelt_topg.f
src/obs_chem.f
src/obs_webmod.f
src/obs_webx.f
src/phreeq_mms.F90
src/potet_hamon_prms.f
src/precip_web.f
src/route_clark.f
src/soltab_prms.f
src/temp_1sta_prms.f
src/top2clark.f
src/topmod_chem.f
src/web_sum.F
src/webmod_res.f
src/xyz_dist_web.f
)

# place [*.F[90]] files into "Source Files"
source_group("Source Files" FILES "../IPhreeqcMMS/IPhreeqc/src/IPhreeqc_interface.F90")
source_group("Source Files" FILES "src/web_sum.F")

# setup include directories
##include_directories("${PROJECT_SOURCE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}/../IPhreeqcMMS/include")
include_directories("${PROJECT_SOURCE_DIR}/../IPhreeqcMMS/IPhreeqc/src")

# executable
add_executable(webmod_f ${webmod_f_SOURCES})
set_target_properties(webmod_f PROPERTIES DEBUG_POSTFIX "d")

# library dependencies
set(EXTRA_LIBS ${EXTRA_LIBS} IPhreeqcMMS)
set(EXTRA_LIBS ${EXTRA_LIBS} IPhreeqc)
set(EXTRA_LIBS ${EXTRA_LIBS} mmf_c)

# link 
target_link_libraries(webmod_f ${EXTRA_LIBS})
