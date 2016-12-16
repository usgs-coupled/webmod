cmake_minimum_required (VERSION 2.8.11)

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

SET(webmod_PROJECTS_ANDREWS
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/Andrews.xlsm
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/webmod.bat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/webmod_gui.bat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/webmod_paramtool.bat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/webmod_print.bat
)

SET(webmod_PROJECTS_ANDREWS_CONTROL
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/control/webmod.control
)

SET(webmod_PROJECTS_ANDREWS_INPUT
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/phreeqc_web_lite.dat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/phreeq_lut
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/webmod.chem.dat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/webmod.hydro.dat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/webmod.params
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/input/webmod.pqi
)

SET(webmod_PROJECTS_ANDREWS_OUTPUT
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/output/webmod.chem.out
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/output/webmod.hydro.out
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/output/webmod.statvar
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews/output/webmod.topout
)

SET(webmod_PROJECTS_ANDREWS_TUTORIAL
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/webmod.bat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/webmod_gui.bat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/webmod_paramtool.bat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/webmod_print.bat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/webmod.sh
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/webmod_gui.sh
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/webmod_paramtool.sh
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/webmod_print.sh
)

SET(webmod_PROJECTS_ANDREWS_TUTORIAL_CONTROL
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/control/webmod.control
)

SET(webmod_PROJECTS_ANDREWS_TUTORIAL_INPUT
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/phreeqc_web_lite.dat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/phreeq_lut
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/webmod.chem.dat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/webmod.hydro.dat
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/webmod.params
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/input/webmod.pqi
)

SET(webmod_PROJECTS_ANDREWS_TUTORIAL_OUTPUT
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/output/webmod.chem.out
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/output/webmod.hydro.out
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/output/webmod.statvar
${CMAKE_SOURCE_DIR}/Distribution/projects/Andrews_tutorial/output/webmod.topout
)

SET(webmod_PROJECTS_DR2
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/dr2.xlsm
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/webmod.bat
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/webmod_gui.bat
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/webmod_paramtool.bat
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/webmod_print.bat
)

SET(webmod_PROJECTS_DR2_CONTROL
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/control/webmod.control
)

SET(webmod_PROJECTS_DR2_INPUT
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/input/phreeqc_web_lite.dat
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/input/phreeq_lut
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/input/webmod.hydro.dat
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/input/webmod.params
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/input/webmod.pqi
)

SET(webmod_PROJECTS_DR2_OUTPUT
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/output/webmod.chem.out
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/output/webmod.hydro.out
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/output/webmod.statvar
${CMAKE_SOURCE_DIR}/Distribution/projects/dr2/output/webmod.topout
)

install(FILES ${webmod_ROOT} DESTINATION .)
install(FILES ${webmod_DOC}  DESTINATION ${CMAKE_INSTALL_DOCDIR})

install(FILES ${webmod_PROJECTS_ANDREWS}         DESTINATION projects/Andrews)
install(FILES ${webmod_PROJECTS_ANDREWS_CONTROL} DESTINATION projects/Andrews/control)
install(FILES ${webmod_PROJECTS_ANDREWS_INPUT}   DESTINATION projects/Andrews/input)
install(FILES ${webmod_PROJECTS_ANDREWS_OUTPUT}  DESTINATION projects/Andrews/output)

install(FILES ${webmod_PROJECTS_ANDREWS_TUTORIAL}         DESTINATION projects/Andrews_tutorial)
install(FILES ${webmod_PROJECTS_ANDREWS_TUTORIAL_CONTROL} DESTINATION projects/Andrews_tutorial/control)
install(FILES ${webmod_PROJECTS_ANDREWS_TUTORIAL_INPUT}   DESTINATION projects/Andrews_tutorial/input)
install(FILES ${webmod_PROJECTS_ANDREWS_TUTORIAL_OUTPUT}  DESTINATION projects/Andrews_tutorial/output)

install(FILES ${webmod_PROJECTS_DR2}         DESTINATION projects/dr2)
install(FILES ${webmod_PROJECTS_DR2_CONTROL} DESTINATION projects/dr2/control)
install(FILES ${webmod_PROJECTS_DR2_INPUT}   DESTINATION projects/dr2/input)
install(FILES ${webmod_PROJECTS_DR2_OUTPUT}  DESTINATION projects/dr2/output)
