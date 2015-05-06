if (WIN32)
  set(WEBMOD_INPUT_ANDCRK_PARAMS         ".\\input\\andcrk.params")
  set(WEBMOD_INPUT_ANDCRK_DAT            ".\\input\\andcrk.dat")
  ###set(WEBMOD_INPUT_PHREEQ_LUT            ".\\input\\phreeq_lut")
  ###set(WEBMOD_INPUT_ANDCRK_DAT_CHEMDAT    ".\\input\\andcrk.dat.chemdat")
  ###set(WEBMOD_INPUT_ANDCRK_DAT_PQI        ".\\input\\andcrk.dat.pqi")
  set(WEBMOD_INPUT_PHREEQC_WEB_LITE_DAT  ".\\input\\phreeqc_web_lite.dat")
  set(WEBMOD_OUTPUT_ANDCRK_OUT           ".\\output\\andcrk.out")
  set(WEBMOD_INPUT_DIR                   ".\\input\\")
  set(WEBMOD_OUTPUT_DIR                  ".\\output\\")
else()
  set(WEBMOD_INPUT_ANDCRK_PARAMS         "./input/andcrk.params")
  set(WEBMOD_INPUT_ANDCRK_DAT            "./input/andcrk.dat")
  ###set(WEBMOD_INPUT_PHREEQ_LUT            "./input/phreeq_lut")
  ###set(WEBMOD_INPUT_ANDCRK_DAT_CHEMDAT    "./input/andcrk.dat.chemdat")
  ###set(WEBMOD_INPUT_ANDCRK_DAT_PQI        "./input/andcrk.dat.pqi")
  set(WEBMOD_INPUT_PHREEQC_WEB_LITE_DAT  "./input/phreeqc_web_lite.dat")
  set(WEBMOD_OUTPUT_ANDCRK_OUT           "./output/andcrk.out")
  set(WEBMOD_INPUT_DIR                   "./input/")
  set(WEBMOD_OUTPUT_DIR                  "./output/")
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
