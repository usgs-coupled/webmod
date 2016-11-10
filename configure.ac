<#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([webmod], [1.0.0-11667], [rmwebb@usgs.gov])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign no-define subdir-objects tar-ustar])
AC_CONFIG_SRCDIR([webmod.d/src/defines.h])

# enable silent rules when available (automake 1.11 or later).
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

AC_CANONICAL_HOST

# IPhreeqc defines
AC_DEFINE([SWIG_SHARED_OBJ], , [Must be set for IPhreeqc])
AC_DEFINE([USE_PHRQ_ALLOC], , [Must be set for IPhreeqc])
AC_DEFINE([USE_GZ], , [Must be set for IPhreeqc])


# Checks for programs.
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_LN_S
AC_PROG_CC
AC_PROG_CPP
AC_PROG_F77
AC_PROG_FC
AC_FC_LIBRARY_LDFLAGS
AC_FC_WRAPPERS
AC_FC_MAIN
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# determine f90 module suffix
AX_F90_MODULE_EXTENSION
F90_MODEXT=mod
if test x$ax_cv_f90_modext != xunknown ; then
  F90_MODEXT=$ax_cv_f90_modext  
fi
AC_SUBST(F90_MODEXT)

# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h memory.h stddef.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memset modf pow sqrt strchr strdup strrchr strstr strtol])

# isfinite
AC_LANG_PUSH(C++)
AC_CACHE_CHECK([for isfinite], ac_cv_isfinite,
 [AC_TRY_LINK([#include <math.h>],
 [double x; int y; y = isfinite(x);],
 ac_cv_isfinite=yes,
 ac_cv_isfinite=no
)])
if test x"$ac_cv_isfinite" = x"yes"; then
  AC_DEFINE(HAVE_ISFINITE, [], [Has isfinite])
fi
AC_LANG_POP(C++)

# finite
AC_LANG_PUSH(C++)
AC_CACHE_CHECK([for finite], ac_cv_finite,
 [AC_TRY_LINK([#include <math.h>],
 [double x; int y; y = finite(x);],
 ac_cv_finite=yes,
 ac_cv_finite=no
)])
if test x"$ac_cv_finite" = x"yes"; then
  AC_DEFINE(HAVE_FINITE, [], [Has finite])
fi
AC_LANG_POP(C++)

# isnan
AC_CHECK_FUNCS([isnan],    , AC_CHECK_LIB(m, isnan))

# set end time
AC_SUBST([WEBMOD_END_YEAR],                  ["1983"])
AC_SUBST([WEBMOD_END_MONTH],                 ["10"])
AC_SUBST([WEBMOD_END_DAY],                   ["6"])

# set paths
AC_SUBST([WEBMOD_INPUT_ANDCRK_PARAMS],         ["./input/andcrk.params"])
AC_SUBST([WEBMOD_INPUT_ANDCRK_DAT],            ["./input/andcrk.dat"])
AC_SUBST([WEBMOD_INPUT_PHREEQC_WEB_LITE_DAT],  ["./input/phreeqc_web_lite.dat"])
AC_SUBST([WEBMOD_INPUT_CHEMDAT_FILE],          ["./input/webmod.chem.dat"])
AC_SUBST([WEBMOD_INPUT_PHREEQ_FILE],           ["./input/webmod.pqi"])
AC_SUBST([WEBMOD_OUTPUT_ANDCRK_OUT],           ["./output/andcrk.out"])
AC_SUBST([WEBMOD_INPUT_DIR],                   ["./input/"])
AC_SUBST([WEBMOD_OUTPUT_DIR],                  ["./output/"])
AC_SUBST([WEBMOD_ANI_OUTPUT],                  ["./output/webmod.aniout"])
AC_SUBST([WEBMOD_STAT_VAR],                    ["./output/andcrk.statvar"])
AC_SUBST([WEBMOD_VAR_INIT],                    ["./input/prms_ic.in"])
AC_SUBST([WEBMOD_SAVE],                        ["./output/andcrk.ics"])  
AC_SUBST([WEBMOD_OPT_OUTPUT],                  ["./output/andcrk.opt"])
AC_SUBST([WEBMOD_PARAM_PRINT],                 ["./output/andcrk.pars"])
AC_SUBST([WEBMOD_SENS_OUTPUT],                 ["./output/andcrk.sens"])
AC_SUBST([WEBMOD_VAR_PRINT],                   ["./output/andcrk.vars"])
AC_SUBST([WEBMOD_GIS_OUTPUT],                  ["./output/andcrk.gis"])
AC_SUBST([WEBMOD_STATS_OUTPUT],                ["./output/andcrk.stats"])
AC_SUBST([WEBMOD_MODEL_CHEMOUT_FILE],          ["./output/webmod.chem.out"])
AC_SUBST([WEBMOD_MODEL_TOPOUT_FILE],           ["./output/webmod.topout"])
AC_SUBST([WEBMOD_EXE_MODEL],                   ["../../bin/webmod_1.0"])

AC_CONFIG_FILES([
  Makefile
  lib/Makefile
  Tests/Makefile
  Tests/Andrews/Makefile
  Tests/Andrews/control/andcrk.control
])

AC_OUTPUT
echo "*************************************"
echo "Summary:"
echo "C++ Compiler                 = $CXX"
echo "CPPFLAGS                     = $CPPFLAGS"
echo "CXXFLAGS                     = $CXXFLAGS"
echo "Fortran Compiler             = $FC"
echo "Fortran Libs                 = $FCLIBS"
echo "LDFLAGS                      = $LDFLAGS"
echo "*************************************"
