# set end_time
set(WEBMOD_END_YEAR                    "1983")
set(WEBMOD_END_MONTH                   "10")
set(WEBMOD_END_DAY                     "6")

# set paths
set(WEBMOD_INPUT_ANDCRK_PARAMS         "./input/andcrk.params")
set(WEBMOD_INPUT_ANDCRK_DAT            "./input/andcrk.dat")
set(WEBMOD_INPUT_PHREEQC_WEB_LITE_DAT  "./input/phreeqc_web_lite.dat")
set(WEBMOD_INPUT_CHEMDAT_FILE          "./input/webmod.chem.dat")
set(WEBMOD_INPUT_PHREEQ_FILE           "./input/webmod.pqi")
set(WEBMOD_OUTPUT_ANDCRK_OUT           "./output/andcrk.out")
set(WEBMOD_INPUT_DIR                   "./input/")
set(WEBMOD_OUTPUT_DIR                  "./output/")
set(WEBMOD_ANI_OUTPUT                  "./output/webmod.aniout")
set(WEBMOD_STAT_VAR                    "./output/andcrk.statvar")
set(WEBMOD_VAR_INIT                    "./input/prms_ic.in")
set(WEBMOD_SAVE                        "./output/andcrk.ics")  
set(WEBMOD_OPT_OUTPUT                  "./output/andcrk.opt")
set(WEBMOD_PARAM_PRINT                 "./output/andcrk.pars")
set(WEBMOD_SENS_OUTPUT                 "./output/andcrk.sens")
set(WEBMOD_VAR_PRINT                   "./output/andcrk.vars")
set(WEBMOD_GIS_OUTPUT                  "./output/andcrk.gis")
set(WEBMOD_STATS_OUTPUT                "./output/andcrk.stats")
set(WEBMOD_MODEL_CHEMOUT_FILE          "./output/webmod.chem.out")
set(WEBMOD_MODEL_TOPOUT_FILE           "./output/webmod.topout")
set(WEBMOD_EXE_MODEL                   "../../bin/webmod_1.0")

# make paths native
file(TO_NATIVE_PATH   ${WEBMOD_INPUT_ANDCRK_PARAMS}         WEBMOD_INPUT_ANDCRK_PARAMS)
file(TO_NATIVE_PATH   ${WEBMOD_INPUT_ANDCRK_DAT}            WEBMOD_INPUT_ANDCRK_DAT)
file(TO_NATIVE_PATH   ${WEBMOD_INPUT_PHREEQC_WEB_LITE_DAT}  WEBMOD_INPUT_PHREEQC_WEB_LITE_DAT)
file(TO_NATIVE_PATH   ${WEBMOD_OUTPUT_ANDCRK_OUT}           WEBMOD_OUTPUT_ANDCRK_OUT)
file(TO_NATIVE_PATH   ${WEBMOD_INPUT_DIR}                   WEBMOD_INPUT_DIR)
file(TO_NATIVE_PATH   ${WEBMOD_OUTPUT_DIR}                  WEBMOD_OUTPUT_DIR)
file(TO_NATIVE_PATH   ${WEBMOD_ANI_OUTPUT}                  WEBMOD_ANI_OUTPUT)
file(TO_NATIVE_PATH   ${WEBMOD_STAT_VAR}                    WEBMOD_STAT_VAR)
file(TO_NATIVE_PATH   ${WEBMOD_VAR_INIT}                    WEBMOD_VAR_INIT)
file(TO_NATIVE_PATH   ${WEBMOD_SAVE}                        WEBMOD_SAVE)
file(TO_NATIVE_PATH   ${WEBMOD_OPT_OUTPUT}                  WEBMOD_OPT_OUTPUT)
file(TO_NATIVE_PATH   ${WEBMOD_PARAM_PRINT}                 WEBMOD_PARAM_PRINT)
file(TO_NATIVE_PATH   ${WEBMOD_SENS_OUTPUT}                 WEBMOD_SENS_OUTPUT)
file(TO_NATIVE_PATH   ${WEBMOD_VAR_PRINT}                   WEBMOD_VAR_PRINT)
file(TO_NATIVE_PATH   ${WEBMOD_GIS_OUTPUT}                  WEBMOD_GIS_OUTPUT)
file(TO_NATIVE_PATH   ${WEBMOD_STATS_OUTPUT}                WEBMOD_STATS_OUTPUT)
file(TO_NATIVE_PATH   ${WEBMOD_EXE_MODEL}                   WEBMOD_EXE_MODEL)

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
  "input/phreeq_lut.in"
  "input/phreeq_lut"
  )
  
# input/andcrk.dat.chemdat  
configure_file(
  "input/andcrk.dat.chemdat"
  "input/andcrk.dat.chemdat"
  )
  
# input/webmod.pqi
configure_file(
  "input/webmod.pqi"
  "input/webmod.pqi"
  )

# input/webmod.chem.dat
configure_file(
  "input/webmod.chem.dat"
  "input/webmod.chem.dat"
  )

# input/phreeqc_web_lite.dat
configure_file(
  "input/phreeqc_web_lite.dat"
  "input/phreeqc_web_lite.dat"
  )


# test andrews creek
add_test(NAME Andrews
         COMMAND webmod
         "-C${PROJECT_BINARY_DIR}/Tests/Andrews/control/andcrk.control"
         )
set_tests_properties(Andrews PROPERTIES PASS_REGULAR_EXPRESSION "Finished")
