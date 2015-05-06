if (WIN32)
  set(WEBMOD_INPUT_ANDCRK_PARAMS ".\\input\\andcrk.params")
else()
  set(WEBMOD_INPUT_ANDCRK_PARAMS "./input/andcrk.params")
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
