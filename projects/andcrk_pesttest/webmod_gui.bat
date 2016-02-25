@ECHO OFF
..\..\bin\webmod -C.\control\webmod.control -print
java -cp ..\..\lib\oui4.jar oui.mms.gui.Mms .\control\webmod.control
ECHO.
ECHO Run complete. Please press enter to continue.
PAUSE>NUL
