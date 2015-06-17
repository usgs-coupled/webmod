echo %time%
@ECHO OFF
..\..\bin\webmod_1.0 -C.\control\webmod.control -print
..\..\bin\webmod_1.0 -C.\control\webmod.control > cmd.out
ECHO.
echo %time%
ECHO Run complete. Please press enter to continue.
PAUSE>NUL