cmake_minimum_required (VERSION 2.8.11)

# use full path on default prefix
if(WEBMOD_USE_ABSOLUTE_PATHS)
  # set jar dir
  set(JARDIR "${CMAKE_INSTALL_FULL_LIBDIR}")

  # set webmod
  set(WEBMOD "${CMAKE_INSTALL_FULL_BINDIR}/webmod")
else()
  # set jar dir
  file(RELATIVE_PATH JARDIR "${CMAKE_INSTALL_FULL_PROJDIR}/Andrews_tutorial" "${CMAKE_INSTALL_FULL_LIBDIR}")

  # set webmod
  file(RELATIVE_PATH WEBMOD "${CMAKE_INSTALL_FULL_PROJDIR}/Andrews_tutorial" "${CMAKE_INSTALL_FULL_BINDIR}/webmod")
endif()

if(WIN32)
  file(TO_NATIVE_PATH ${JARDIR} JARDIR)
  file(TO_NATIVE_PATH ${WEBMOD} WEBMOD)
  if(WEBMOD_USE_ABSOLUTE_PATHS)
    set(JARDIR "\"${JARDIR}\"")
    set(WEBMOD "\"${WEBMOD}\"")
  endif()
endif()

foreach(dir
    Andrews
    Andrews_tutorial
    dr2
    )
  # webmod.bat
  configure_file(
    "${CMAKE_SOURCE_DIR}/Distribution/projects/${dir}/webmod.bat.in"
    "${CMAKE_BINARY_DIR}/Distribution/projects/${dir}/webmod.bat"
    )
  # webmod_gui.bat
  configure_file(
    "${CMAKE_SOURCE_DIR}/Distribution/projects/${dir}/webmod_gui.bat.in"
    "${CMAKE_BINARY_DIR}/Distribution/projects/${dir}/webmod_gui.bat"
    )
  # webmod_paramtool.bat
  configure_file(
    "${CMAKE_SOURCE_DIR}/Distribution/projects/${dir}/webmod_paramtool.bat.in"
    "${CMAKE_BINARY_DIR}/Distribution/projects/${dir}/webmod_paramtool.bat"
    )
  # webmod_print.bat
  configure_file(
    "${CMAKE_SOURCE_DIR}/Distribution/projects/${dir}/webmod_print.bat.in"
    "${CMAKE_BINARY_DIR}/Distribution/projects/${dir}/webmod_print.bat"
    )

    
  # webmod.sh
  configure_file(
    "${CMAKE_SOURCE_DIR}/Distribution/projects/${dir}/webmod.sh.in"
    "${CMAKE_BINARY_DIR}/Distribution/projects/${dir}/webmod.sh"
    )
  # webmod_gui.sh
  configure_file(
    "${CMAKE_SOURCE_DIR}/Distribution/projects/${dir}/webmod_gui.sh.in"
    "${CMAKE_BINARY_DIR}/Distribution/projects/${dir}/webmod_gui.sh"
    )
  # webmod_paramtool.sh
  configure_file(
    "${CMAKE_SOURCE_DIR}/Distribution/projects/${dir}/webmod_paramtool.sh.in"
    "${CMAKE_BINARY_DIR}/Distribution/projects/${dir}/webmod_paramtool.sh"
    )
  # webmod_print.sh
  configure_file(
    "${CMAKE_SOURCE_DIR}/Distribution/projects/${dir}/webmod_print.sh.in"
    "${CMAKE_BINARY_DIR}/Distribution/projects/${dir}/webmod_print.sh"
    )


  # webmod.control
  configure_file(
    "${CMAKE_SOURCE_DIR}/Distribution/projects/${dir}/control/webmod.control.in"
    "${CMAKE_BINARY_DIR}/Distribution/projects/${dir}/control/webmod.control"
    )
endforeach()

SET(webmod_ROOT
  ${CMAKE_SOURCE_DIR}/README_WEBMOD.TXT
  ${CMAKE_SOURCE_DIR}/Distribution/README.txt
  )

SET(webmod_DOC
  ${CMAKE_SOURCE_DIR}/Distribution/doc/Webb_WEBMOD_Body_2016Aug19.docx
  ${CMAKE_SOURCE_DIR}/Distribution/doc/Webb_WEBMOD_Tables_2016August19.docm
  ${CMAKE_SOURCE_DIR}/Distribution/doc/WEBMOD_figs_1_12_2016August19.pdf
  ${CMAKE_SOURCE_DIR}/Distribution/doc/WEBMOD_figs_13_52_2016August19.pdf
  )

if (WIN32)
  SET(webmod_PROJECTS_ANDREWS
    "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/Andrews.xlsm"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews/webmod.bat"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews/webmod_gui.bat"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews/webmod_paramtool.bat"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews/webmod_print.bat"
    )
else()
  SET(webmod_PROJECTS_ANDREWS
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews/webmod.sh"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews/webmod_gui.sh"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews/webmod_paramtool.sh"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews/webmod_print.sh"
    )
endif()

SET(webmod_PROJECTS_ANDREWS_CONTROL
  "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews/control/webmod.control"
)

SET(webmod_PROJECTS_ANDREWS_INPUT
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/phreeqc_web_lite.dat"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/phreeq_lut"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/webmod.chem.dat"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/webmod.hydro.dat"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/webmod.params"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/webmod.pqi"
)

SET(webmod_PROJECTS_ANDREWS_OUTPUT
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/output/webmod.chem.out"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/output/webmod.hydro.out"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/output/webmod.statvar"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/output/webmod.topout"
)

if (WIN32)
  SET(webmod_PROJECTS_ANDREWS_TUTORIAL
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews_tutorial/webmod.bat"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews_tutorial/webmod_gui.bat"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews_tutorial/webmod_paramtool.bat"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews_tutorial/webmod_print.bat"
  )
else()
  SET(webmod_PROJECTS_ANDREWS_TUTORIAL
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews_tutorial/webmod.sh"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews_tutorial/webmod_gui.sh"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews_tutorial/webmod_paramtool.sh"
    "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews_tutorial/webmod_print.sh"
  )
endif()

SET(webmod_PROJECTS_ANDREWS_TUTORIAL_CONTROL
  "${CMAKE_BINARY_DIR}/Distribution/projects/Andrews_tutorial/control/webmod.control"
  )

SET(webmod_PROJECTS_ANDREWS_TUTORIAL_INPUT
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/phreeqc_web_lite.dat"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/phreeq_lut"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/webmod.chem.dat"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/webmod.hydro.dat"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/webmod.params"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/webmod.pqi"
  )

SET(webmod_PROJECTS_ANDREWS_TUTORIAL_OUTPUT
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/output/webmod.chem.out"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/output/webmod.hydro.out"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/output/webmod.statvar"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/output/webmod.topout"
  )

if (WIN32)
  SET(webmod_PROJECTS_DR2
    "${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/dr2.xlsm"
    "${CMAKE_BINARY_DIR}/Distribution/projects/dr2/webmod.bat"
    "${CMAKE_BINARY_DIR}/Distribution/projects/dr2/webmod_gui.bat"
    "${CMAKE_BINARY_DIR}/Distribution/projects/dr2/webmod_paramtool.bat"
    "${CMAKE_BINARY_DIR}/Distribution/projects/dr2/webmod_print.bat"
    )
else()
  SET(webmod_PROJECTS_DR2
    "${CMAKE_BINARY_DIR}/Distribution/projects/dr2/webmod.sh"
    "${CMAKE_BINARY_DIR}/Distribution/projects/dr2/webmod_gui.sh"
    "${CMAKE_BINARY_DIR}/Distribution/projects/dr2/webmod_paramtool.sh"
    "${CMAKE_BINARY_DIR}/Distribution/projects/dr2/webmod_print.sh"
    )
endif()

SET(webmod_PROJECTS_DR2_CONTROL
  "${CMAKE_BINARY_DIR}/Distribution/projects/dr2/control/webmod.control"
  )

SET(webmod_PROJECTS_DR2_INPUT
  "${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/input/phreeqc_web_lite.dat"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/input/phreeq_lut"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/input/webmod.hydro.dat"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/input/webmod.params"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/input/webmod.pqi"
  )

SET(webmod_PROJECTS_DR2_OUTPUT
  "${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/output/webmod.chem.out"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/output/webmod.hydro.out"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/output/webmod.statvar"
  "${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/output/webmod.topout"
)

if (NOT WEBMOD_USE_ABSOLUTE_PATHS OR WIN32)
  install(FILES ${webmod_ROOT} DESTINATION ".")
else()
  install(FILES ${webmod_ROOT} DESTINATION "${CMAKE_INSTALL_DOCDIR}")
endif()
install(FILES ${webmod_DOC} DESTINATION "${CMAKE_INSTALL_DOCDIR}")

if (WIN32)
  install(FILES ${webmod_PROJECTS_ANDREWS} DESTINATION "${CMAKE_INSTALL_PROJDIR}/Andrews")
else()
  install(PROGRAMS ${webmod_PROJECTS_ANDREWS} DESTINATION "${CMAKE_INSTALL_PROJDIR}/Andrews")
  install(FILES "${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/Andrews.xlsm" DESTINATION "${CMAKE_INSTALL_PROJDIR}/Andrews")
endif()
install(FILES ${webmod_PROJECTS_ANDREWS_CONTROL} DESTINATION "${CMAKE_INSTALL_PROJDIR}/Andrews/control")
install(FILES ${webmod_PROJECTS_ANDREWS_INPUT}   DESTINATION "${CMAKE_INSTALL_PROJDIR}/Andrews/input")
install(FILES ${webmod_PROJECTS_ANDREWS_OUTPUT}  DESTINATION "${CMAKE_INSTALL_PROJDIR}/Andrews/output")

if (WIN32)
  install(FILES ${webmod_PROJECTS_ANDREWS_TUTORIAL} DESTINATION "${CMAKE_INSTALL_PROJDIR}/Andrews_tutorial")
else()
  install(PROGRAMS ${webmod_PROJECTS_ANDREWS_TUTORIAL} DESTINATION "${CMAKE_INSTALL_PROJDIR}/Andrews_tutorial")
endif()
install(FILES ${webmod_PROJECTS_ANDREWS_TUTORIAL_CONTROL} DESTINATION "${CMAKE_INSTALL_PROJDIR}/Andrews_tutorial/control")
install(FILES ${webmod_PROJECTS_ANDREWS_TUTORIAL_INPUT}   DESTINATION "${CMAKE_INSTALL_PROJDIR}/Andrews_tutorial/input")
install(FILES ${webmod_PROJECTS_ANDREWS_TUTORIAL_OUTPUT}  DESTINATION "${CMAKE_INSTALL_PROJDIR}/Andrews_tutorial/output")

if (WIN32)
  install(FILES ${webmod_PROJECTS_DR2} DESTINATION "${CMAKE_INSTALL_PROJDIR}/dr2")
else()
  install(PROGRAMS ${webmod_PROJECTS_DR2} DESTINATION "${CMAKE_INSTALL_PROJDIR}/dr2")
  install(FILES "${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/dr2.xlsm" DESTINATION "${CMAKE_INSTALL_PROJDIR}/dr2")
endif()
install(FILES ${webmod_PROJECTS_DR2_CONTROL} DESTINATION "${CMAKE_INSTALL_PROJDIR}/dr2/control")
install(FILES ${webmod_PROJECTS_DR2_INPUT}   DESTINATION "${CMAKE_INSTALL_PROJDIR}/dr2/input")
install(FILES ${webmod_PROJECTS_DR2_OUTPUT}  DESTINATION "${CMAKE_INSTALL_PROJDIR}/dr2/output")
