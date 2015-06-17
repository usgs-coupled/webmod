@echo off
REM Sychronize latest WEBMOD executable, create the tsproc and PEST command files, and start beopest as master
REM Help printed if command entered with no arguments
IF %1.==. GOTO usage
IF %2.==. GOTO NoCase

set HOME=%cd%
set HOME_SED=%HOME:\=/%

REM Set variables
set nodes=%1
set pst=%2
set BASIN=andcrk

set "INPUT_DIR=%HOME%\input"
set INPUT_DIR_SED=%INPUT_DIR:\=/%
set "PEST_FILES_DIR=%HOME%\pest_files"
set PEST_FILES_DIR_SED=%PEST_FILES_DIR:\=/%
set "CONTROL_DIR=%HOME%\control"
set CONTROL_DIR_SED=%CONTROL:\=/%
set "PROJECT_DIR=%HOME%\pest_run_dir"
set PROJECT_DIR_SED=%PROJECT_DIR:\=/%
set "PEST_BIN_DIR=%HOME%\..\..\bin"
set PEST_BIN_DIR_SED=%PEST_BIN_DIR:\=/%
set "TSPROC_BIN_DIR=%HOME%\..\..\bin"
set TSPROC_BIN_DIR_SED=%TSPROC_BIN_DIR:\=/%

REM set port number
set PORT=4004 

REM setup working directory PROJECT_DIR=pest_run_dir
if exist %PROJECT_DIR% rmdir /s/q %PROJECT_DIR%
mkdir %PROJECT_DIR%

REM Sed to make basin_tsproc.dat 
call %HOME%\sed.bat @PROJECT_DIR@ %PROJECT_DIR_SED%/ %PEST_FILES_DIR%\%BASIN%_tsproc.tpl > %PROJECT_DIR%\%BASIN%_tsproc.dat

REM Sed to make pest_webmod.bat in PROJECT_DIR
call %HOME%\sed.bat @PROJECT_DIR@ %PROJECT_DIR_SED%/   %PEST_FILES_DIR%\pest_webmod.bat.tpl > %PROJECT_DIR%\pest_webmod.bat1
call %HOME%\sed.bat @PEST_BIN_DIR@ %PEST_BIN_DIR_SED%/ %PROJECT_DIR%\pest_webmod.bat1 >       %PROJECT_DIR%\pest_webmod.bat2
call %HOME%\sed.bat @TSPROC_BIN_DIR@ %TSPROC_BIN_DIR%/ %PROJECT_DIR%\pest_webmod.bat2 >       %PROJECT_DIR%\pest_webmod.bat3
call %HOME%\sed.bat @DEL@ DEL                          %PROJECT_DIR%\pest_webmod.bat3 >       %PROJECT_DIR%\pest_webmod.bat
DEL %PROJECT_DIR%\pest_webmod.bat1 %PROJECT_DIR%\pest_webmod.bat2 %PROJECT_DIR%\pest_webmod.bat3

REM Copy files to tsproc directory
copy %PEST_FILES_DIR%\*.ssf            	        %PROJECT_DIR%
copy %PEST_FILES_DIR%\pest_groups*.txt 	        %PROJECT_DIR%
copy %PEST_FILES_DIR%\pest_params*.txt 	        %PROJECT_DIR%
copy %PEST_FILES_DIR%\tsproc.in         	%PROJECT_DIR%
copy %PEST_FILES_DIR%\%BASIN%.statvar.full	%PROJECT_DIR%\%BASIN%.statvar
copy %PEST_FILES_DIR%\par2par_%BASIN%.tpl	%PROJECT_DIR%

REM run tsproc
cd %PROJECT_DIR%
%PEST_BIN_DIR%\tsproc.exe < .\tsproc.in

REM move output from tsproc to working directory
COPY %TSPROC_BIN_DIR%\webmod_1.0.exe  %PROJECT_DIR%

REM edit the control file
call %HOME%\sed.bat ./input . %CONTROL_DIR%\%BASIN%.control >    %PROJECT_DIR%\%BASIN%.control1
call %HOME%\sed.bat ./output . %PROJECT_DIR%\%BASIN%.control1 >  %PROJECT_DIR%\%BASIN%.control2
call %HOME%\sed.bat ../../bin . %PROJECT_DIR%\%BASIN%.control2 >  %PROJECT_DIR%\%BASIN%.control
DEL %PROJECT_DIR%\%BASIN%.control1 %PROJECT_DIR%\%BASIN%.control2


REM change tsproc mode from "pest_prep" to "model_run"
call %HOME%\sed.bat "CONTEXT pest_prep" "CONTEXT model_run" %PROJECT_DIR%\%BASIN%_tsproc.dat > %PROJECT_DIR%\xxx
MOVE  %PROJECT_DIR%\xxx %PROJECT_DIR%\%BASIN%_tsproc.dat


REM Line 6 of the pst file RLAMBDA1 .....	
rem sed -i "6s/.*/10.0  -3.0    0.3    0.03     -10  999  LAMFORGIVE/" %PROJECT_DIR_PATH%\%BASIN%.pst
REM Line 7 of the pst file RELPARMAX FACPARMAX FACORIG.....	
rem sed -i "7s/.*/0.2   2.0   1.0e-3/" %PROJECT_DIR_PATH%\%BASIN%.pst
REM Line 9 of the pst file. NOPTMAX Max # of optimizations. Default is 30, set to 0 for single run with phi contributions, 1 for sensitivities, or a small number to test PEST loops.
rem sed -i "9s/.*/30   .005  4   4  .005   4/" %PROJECT_DIR_PATH%\%BASIN%.pst
REM Line 13. SVD block MAXSING EIGTHRESH. Replace MAXSING the maximum number of adjustable variables (number of singlular valuess at which truncation occurs)
rem sed -i "13s/.*/19 5e-7/" %PROJECT_DIR_PATH%\%BASIN%.pst
REM Line 14 EIGWRITE. 0 if not using SVD output file
rem sed -i "14s/.*/1/" %PROJECT_DIR_PATH%\%BASIN%.pst
REM remove sed files
rem rm sed*

REM Edit .pst file for parameter estimation
SETLOCAL ENABLEEXTENSIONS
SETLOCAL ENABLEDELAYEDEXPANSION
set /a _COUNTER=0
set pstfile=%PROJECT_DIR%\%BASIN%.pst
set tmppest=%PROJECT_DIR%\tmp.pst
for /f "tokens=1,* delims=]" %%A in ('"type %pstfile%|find /n /v """') do (
    set /a "COUNTER+=1"
    set "line=%%B"  
    if "!COUNTER!" == "1" (
        ECHO !line! > %tmppest%
    ) else if "!COUNTER!" == "6" (
        ECHO 10.0  -3.0    0.3    0.03     -10  999  LAMFORGIVE >> %tmppest%
    ) else if "!COUNTER!" == "7" ( 
        ECHO 0.2   2.0   1.0e-3 >> %tmppest%
    ) else if "!COUNTER!" == "9" ( 
        ECHO 30   .005  4   4  .005   4 >> %tmppest%
    ) else if "!COUNTER!" == "13" ( 
        ECHO 19 5e-7 >> %tmppest%
    ) else if "!COUNTER!" == "14" ( 
        ECHO 1 >> %tmppest%
    ) else (
        ECHO !line! >> %tmppest%
    )
)
MOVE %tmppest% %pstfile%

REM Run parallel pest Master
cd %PROJECT_DIR%
start "Master" cmd /k "call %PEST_BIN_DIR%\beopest64.exe %PROJECT_DIR%\%pst% /H :%PORT% & cd .."

REM wait 2 seconds for master to come up before starting workers
sleep 2

REM make tmp dirs and run parallel pest workers
set MASTER=%COMPUTERNAME%
for /l %%X in (1, 1, %nodes%) do (
	if exist %PROJECT_DIR%\tmpest%%X rmdir /s/q %PROJECT_DIR%\tmpest%%X 
	mkdir %PROJECT_DIR%\tmpest%%X
	cd %PROJECT_DIR%\tmpest%%X
	copy %PEST_FILES_DIR%\params_%BASIN%.tpl .\
	copy %PEST_FILES_DIR%\pqi_%BASIN%.tpl .\
	copy %INPUT_DIR%\%BASIN%.dat .\
	copy %INPUT_DIR%\%BASIN%.dat.chemdat .\
	copy %INPUT_DIR%\phreeq_lut .\
	copy %INPUT_DIR%\phreeqc_web_lite.dat .\
	copy %PROJECT_DIR%\pest_webmod.bat .\
	copy %PROJECT_DIR%\%BASIN%.ins .\
	copy %PROJECT_DIR%\%BASIN%_tsproc.dat .\
	START /B %PEST_BIN_DIR%\beopest64.exe %PROJECT_DIR%\%pst% /H %MASTER%:%PORT%
	cd %PROJECT_DIR%
)
cd %PROJECT_DIR_PATH%\..

REM Workers have been started.
GOTO :EOF

:usage
echo.
echo usage: runpest n case.pst
echo where  n is number of slaves to start, and
echo        case.pst is the name of the PEST control file.
GOTO :EOF

:NoCase
echo.
echo Missing name of PEST control file
echo.
echo usage: runslave n case.pst
echo where  n is number of slaves to start, and
echo        case.pst is the name of the PEST control file.
:EOF
cd %PROJECT_DIR_PATH%\..


