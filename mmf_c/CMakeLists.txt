cmake_minimum_required (VERSION 2.8)

project(mmf_c C)

# fortran name mangling
include(FortranCInterface)
FortranCInterface_HEADER(FC.h MACRO_NAMESPACE "FC_" SYMBOLS
control_string_
dattim
decldim
declfix
declmodule
declparam
declpri
declvar
decl_control
deltim
djulian
dpdble
dpint4_
dpreal_
dpstr
getdataname
getdim
getdimname
getoutname
getparam
getstep
getvar
isleap
julian
opstr
putvar
readvar
setdims_
##setdims
##call_modules_
call_modules
)

# set "Character Set" to Unicode
# this may be unnecessary
add_definitions(-D_UNICODE)
add_definitions(-DUNICODE)

# turn off "The POSIX name for this item is deprecated" warning
add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)

SET(mmf_c_SOURCES
src/alloc_space.c
src/batch_run.c
src/batch_run_functions.c
src/build_lists.c
src/call_modules.c
src/call_setdims.c
src/check_vars.c
src/control_addr.c
src/control_array.c
src/control_var.c
src/create_vstats.c
src/decl_control.c
src/decldim.c
src/declparam.c
src/declvar.c
src/defs.h
src/dim_addr.c
src/dprint.c
src/free_vstats.c
src/get_elem_add.c
src/get_times.c
src/getdim.c
src/getdimname.c
src/getparam.c
src/getvar.c
src/graph_single_run.c
src/julconvert.c
src/julday.c
src/load_param.c
src/mmf.c
src/oprint.c
src/param_addr.c
src/parse_args.c
src/print_model_info.c
src/print_params.c
src/print_vars.c
src/putvar.c
src/read_control.c
src/read_datainfo.c
src/read_line.c
src/read_params.c
src/read_vars.c
src/readvar.c
src/reset_dim.c
src/save_params.c
src/save_vars.c
src/setup_cont.c
src/sort_dims.c
src/sort_params.c
src/sort_vars.c
src/stats.c
src/str_to_vals.c
src/timing.c
src/umalloc_etc.c
src/uprint.c
src/var_addr.c
src/write_vstats.c
)


# setup include directories
include_directories("${PROJECT_BINARY_DIR}")
##include_directories("${PROJECT_SOURCE_DIR}/IPhreeqc/src")
##include_directories("${PROJECT_SOURCE_DIR}/IPhreeqc/src/phreeqcpp")
##include_directories("${PROJECT_SOURCE_DIR}/IPhreeqc/src/phreeqcpp/PhreeqcKeywords")

# library
add_library(mmf_c ${LIB_TYPE} ${mmf_c_SOURCES})
set_target_properties(mmf_c PROPERTIES DEBUG_POSTFIX "d")
