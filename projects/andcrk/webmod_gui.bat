@ECHO OFF
..\..\bin\webmod_1.0 -C.\control\andcrk.control -print
java -cp ..\..\lib\oui3.jar;..\..\lib\AbsoluteLayout.jar;..\..\lib\jcommon-1.0.12.jar;..\..\lib\jfreechart-1.0.9.jar oui.mms.gui.Mms .\control\webmod.control
ECHO.
ECHO Run complete. Please press enter to continue.
PAUSE>NUL
