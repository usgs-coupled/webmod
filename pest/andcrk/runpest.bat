@echo off
REM Sychronize latest WEBMOD executable, create the tsproc and PEST command files, and start beopest as master
REM Help printed if command entered with no arguments
IF %1.==. GOTO :USAGE

REM Set variables
set NODES=%1
set PST=webmod.pst
set HOME=%cd%
set HOME_SED=%HOME:\=/%
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

REM Sed to make tsproc.dat 
call %HOME%\sed.bat @PROJECT_DIR@ %PROJECT_DIR_SED%/ %PEST_FILES_DIR%\tsproc.dat.tpl > %PROJECT_DIR%\tsproc.dat

REM Sed to make pest_webmod.bat in PROJECT_DIR
call %HOME%\sed.bat @PROJECT_DIR@ %PROJECT_DIR_SED%/       %PEST_FILES_DIR%\pest_webmod.bat.tpl > %PROJECT_DIR%\pest_webmod.bat1
call %HOME%\sed.bat @PEST_BIN_DIR@ %PEST_BIN_DIR_SED%/     %PROJECT_DIR%\pest_webmod.bat1 >       %PROJECT_DIR%\pest_webmod.bat2
call %HOME%\sed.bat @TSPROC_BIN_DIR@ %TSPROC_BIN_DIR_SED%/ %PROJECT_DIR%\pest_webmod.bat2 >       %PROJECT_DIR%\pest_webmod.bat3
call %HOME%\sed.bat @DEL@ DEL                              %PROJECT_DIR%\pest_webmod.bat3 >       %PROJECT_DIR%\pest_webmod.bat
DEL %PROJECT_DIR%\pest_webmod.bat1 
DEL %PROJECT_DIR%\pest_webmod.bat2 
DEL %PROJECT_DIR%\pest_webmod.bat3

REM Copy files to tsproc directory
copy %PEST_FILES_DIR%\*.ssf            	        %PROJECT_DIR%
copy %PEST_FILES_DIR%\pest_groups*.txt 	        %PROJECT_DIR%
copy %PEST_FILES_DIR%\pest_params*.txt 	        %PROJECT_DIR%
copy %PEST_FILES_DIR%\tsproc.in         	%PROJECT_DIR%
copy %PEST_FILES_DIR%\webmod.statvar	        %PROJECT_DIR%\webmod.statvar
copy %PEST_FILES_DIR%\par2par.dat.tpl	        %PROJECT_DIR%

REM run tsproc
cd %PROJECT_DIR%
%PEST_BIN_DIR%\tsproc.exe < .\tsproc.in

REM move output from tsproc to working directory
COPY %TSPROC_BIN_DIR%\webmod.exe  %PROJECT_DIR%

REM edit the control file
call %HOME%\sed.bat   ./input  .  %CONTROL_DIR%\webmod.control  >  %PROJECT_DIR%\webmod.control1
call %HOME%\sed.bat  ./output  .  %PROJECT_DIR%\webmod.control1 >  %PROJECT_DIR%\webmod.control2
call %HOME%\sed.bat ../../bin  .  %PROJECT_DIR%\webmod.control2 >  %PROJECT_DIR%\webmod.control
DEL %PROJECT_DIR%\webmod.control1  
DEL %PROJECT_DIR%\webmod.control2

REM change tsproc mode from "pest_prep" to "model_run"
call %HOME%\sed.bat "CONTEXT pest_prep" "CONTEXT model_run" %PROJECT_DIR%\tsproc.dat > %PROJECT_DIR%\xxx
MOVE  %PROJECT_DIR%\xxx %PROJECT_DIR%\tsproc.dat

REM Edit .pst file for parameter estimation
SETLOCAL ENABLEEXTENSIONS
SETLOCAL ENABLEDELAYEDEXPANSION
set /a COUNTER=0
set pstfile=%PROJECT_DIR%\webmod.pst
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

REM make tmp dirs and run parallel pest workers
set MASTER=%COMPUTERNAME%
for /l %%X in (1, 1, %NODES%) do (
	if exist %PROJECT_DIR%\tmpest%%X rmdir /s/q %PROJECT_DIR%\tmpest%%X 
	mkdir %PROJECT_DIR%\tmpest%%X
	cd %PROJECT_DIR%\tmpest%%X
	copy %PEST_FILES_DIR%\webmod.params.tpl .\
	copy %PEST_FILES_DIR%\webmod.pqi.tpl    .\
	copy %INPUT_DIR%\webmod.hydro.dat       .\
	copy %INPUT_DIR%\webmod.chem.dat        .\
	copy %INPUT_DIR%\phreeq_lut             .\
	copy %INPUT_DIR%\phreeqc_web_lite.dat   .\
	copy %PROJECT_DIR%\pest_webmod.bat      .\
	copy %PROJECT_DIR%\tsproc.sim.out.ins   .\
	copy %PROJECT_DIR%\tsproc.dat           .\
	START /B %PEST_BIN_DIR%\beopest64.exe %PROJECT_DIR%\%PST% /H %MASTER%:%PORT% & cd %PROJECT_DIR%\..
	cd %PROJECT_DIR%\..
)

REM /L requires a run in PROJECT_DIR
cd %PROJECT_DIR%
copy %PEST_FILES_DIR%\webmod.params.tpl .\
copy %PEST_FILES_DIR%\webmod.pqi.tpl    .\
copy %INPUT_DIR%\webmod.hydro.dat       .\
copy %INPUT_DIR%\webmod.chem.dat        .\
copy %INPUT_DIR%\phreeq_lut             .\
copy %INPUT_DIR%\phreeqc_web_lite.dat   .\

REM Run parallel pest Master
cd %PROJECT_DIR%
%PEST_BIN_DIR%\beopest64.exe %PROJECT_DIR%\%PST% /H /L :%PORT% 

REM tidy up
IF "%ERRORLEVEL%" == "0" (
    if exist %PROJECT_DIR%\..\pest_results rmdir /s/q %PROJECT_DIR%\..\pest_results
    mkdir %PROJECT_DIR%\..\pest_results
    cd %PROJECT_DIR%
    copy tsproc.sim.out.ins ..\pest_results
    copy webmod.jac         ..\pest_results
    copy webmod.jco         ..\pest_results
    copy webmod.jst         ..\pest_results
    copy webmod.mtt         ..\pest_results
    copy webmod.par         ..\pest_results
    copy webmod.prf         ..\pest_results
    copy webmod.pst         ..\pest_results
    copy webmod.rec         ..\pest_results
    copy webmod.rei         ..\pest_results
    copy webmod.res         ..\pest_results
    copy webmod.rmr         ..\pest_results
    copy webmod.rsd         ..\pest_results
    copy webmod.rst         ..\pest_results
    copy webmod.sen         ..\pest_results
    copy webmod.seo         ..\pest_results
    copy webmod.svd         ..\pest_results
    copy webmod.pqi         ..\pest_results
    copy webmod.params      ..\pest_results
    copy webmod.pqi         ..\input\webmod.pqi.pest
    copy webmod.params      ..\input\webmod.params.pest
    copy tsproc.dat         ..\pest_results
    cd   %PROJECT_DIR%\..
REM
REM run sensitivity plots here
REM
) else (
    echo Beopest failed.
)

REM Done.
GOTO :FINISH_UP

:USAGE
echo.
echo usage: runpest n 
echo where  n is number of workers to start
GOTO :EOF

:FINISH_UP

cd %PROJECT_DIR%\..


