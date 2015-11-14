echo %time%
@ECHO OFF
..\..\bin\webmod.exe -C.\control\webmod.control -print
..\..\bin\webmod.exe -C.\control\webmod.control > webmod.log
ECHO.
echo %time%
ECHO Run complete. Please press enter to continue.
PAUSE>NUL