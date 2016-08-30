cmake_minimum_required (VERSION 2.8.11)

SET(webmod_LIB
oui4.jar
README.TXT
)

SET(webmod_LIB_LIB
lib/jfreechart-1.0.13.jar
lib/jhall.jar
lib/swing-layout-1.0.4.jar
lib/AbsoluteLayout.jar
lib/csvutils_old.jar
lib/jcommon-1.0.16.jar
)

install (FILES ${webmod_LIB} DESTINATION ${CMAKE_INSTALL_LIBDIR})
install (FILES ${webmod_LIB_LIB} DESTINATION ${CMAKE_INSTALL_LIBDIR}/lib)
