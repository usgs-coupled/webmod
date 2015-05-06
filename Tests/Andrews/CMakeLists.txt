set(WEBMOD_END_YEAR                    "1983")
set(WEBMOD_END_MONTH                   "10")
set(WEBMOD_END_DAY                     "6")

if (WIN32)
  set(WEBMOD_INPUT_ANDCRK_PARAMS         ".\\input\\andcrk.params")
  set(WEBMOD_INPUT_ANDCRK_DAT            ".\\input\\andcrk.dat")
  set(WEBMOD_INPUT_PHREEQC_WEB_LITE_DAT  ".\\input\\phreeqc_web_lite.dat")
  set(WEBMOD_OUTPUT_ANDCRK_OUT           ".\\output\\andcrk.out")
  set(WEBMOD_INPUT_DIR                   ".\\input\\")
  set(WEBMOD_OUTPUT_DIR                  ".\\output\\")
  set(WEBMOD_ANI_OUTPUT                  ".\\output\\animation.out")
  set(WEBMOD_STAT_VAR                    ".\\output\\andcrk.statvar")
  set(WEBMOD_VAR_INIT                    ".\\input\\prms_ic.in")
  set(WEBMOD_SAVE                        ".\\output\\andcrk.ics")
  set(WEBMOD_OPT_OUTPUT                  ".\\output\\andcrk.opt")
  set(WEBMOD_PARAM_PRINT                 ".\\output\\andcrk.pars")
  set(WEBMOD_SENS_OUTPUT                 ".\\output\\andcrk.sens")
  set(WEBMOD_VAR_PRINT                   ".\\output\\andcrk.vars")
  set(WEBMOD_GIS_OUTPUT                  ".\\output\\andcrk.gis")
  set(WEBMOD_STATS_OUTPUT                ".\\output\\andcrk.stats")
  set(WEBMOD_EXE_MODEL                   "..\\..\\bin\\webmod_1.0")
else()
  set(WEBMOD_INPUT_ANDCRK_PARAMS         "./input/andcrk.params")
  set(WEBMOD_INPUT_ANDCRK_DAT            "./input/andcrk.dat")
  set(WEBMOD_INPUT_PHREEQC_WEB_LITE_DAT  "./input/phreeqc_web_lite.dat")
  set(WEBMOD_OUTPUT_ANDCRK_OUT           "./output/andcrk.out")
  set(WEBMOD_INPUT_DIR                   "./input/")
  set(WEBMOD_OUTPUT_DIR                  "./output/")
  set(WEBMOD_ANI_OUTPUT                  "./output/animation.out")
  set(WEBMOD_STAT_VAR                    "./output/andcrk.statvar")
  set(WEBMOD_VAR_INIT                    "./input/prms_ic.in")
  set(WEBMOD_SAVE                        "./output/andcrk.ics")  
  set(WEBMOD_OPT_OUTPUT                  "./output/andcrk.opt")
  set(WEBMOD_PARAM_PRINT                 "./output/andcrk.pars")
  set(WEBMOD_SENS_OUTPUT                 "./output/andcrk.sens")
  set(WEBMOD_VAR_PRINT                   "./output/andcrk.vars")
  set(WEBMOD_GIS_OUTPUT                  "./output/andcrk.gis")
  set(WEBMOD_STATS_OUTPUT                "./output/andcrk.stats")
  set(WEBMOD_EXE_MODEL                   "../../bin/webmod_1.0")
endif()

# create output directory
file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/Tests/Andrews/output")

# control/andcrk.control
configure_file(
  "control/andcrk.control.in"
  "control/andcrk.control"
  )

# input/andcrk.params
configure_file(
  "input/andcrk.params"
  "input/andcrk.params"
  )
  
# input/andcrk.dat
configure_file(
  "input/andcrk.dat"
  "input/andcrk.dat"
  )
  
# input/phreeq_lut
configure_file(
  "input/phreeq_lut"
  "input/phreeq_lut"
  )
  
# input/andcrk.dat.chemdat  
configure_file(
  "input/andcrk.dat.chemdat"
  "input/andcrk.dat.chemdat"
  )
  
# input/andcrk.dat.pqi  
configure_file(
  "input/andcrk.dat.pqi"
  "input/andcrk.dat.pqi"
  )
  
# input/phreeqc_web_lite.dat
configure_file(
  "input/phreeqc_web_lite.dat"
  "input/phreeqc_web_lite.dat"
  )


# test run
add_test(NAME Andrews
         COMMAND webmod_f
         "-C${PROJECT_BINARY_DIR}/Tests/Andrews/control/andcrk.control"
         )
set_tests_properties(Andrews PROPERTIES PASS_REGULAR_EXPRESSION "Finished")         
