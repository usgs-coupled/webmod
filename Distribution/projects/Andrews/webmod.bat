echo %time%
@ECHO OFF
..\..\bin\webmod -C.\control\webmod.control -print
..\..\bin\webmod -C.\control\webmod.control
ECHO.
echo %time%
ECHO Run complete. Please press enter to continue.
PAUSE>NUL