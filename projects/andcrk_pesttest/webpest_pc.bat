@echo off
SETLOCAL

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
set "OUTPUT_DIR=%HOME%\output"
set OUTPUT_DIR_SED=%OUTPUT_DIR:\=/%
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

REM check hydro files
set errors=0
for %%i in ( 
    %INPUT_DIR%\webmod.hydro.dat
    %INPUT_DIR%\webmod.params
    %CONTROL_DIR%\webmod.control
    %PEST_FILES_DIR%\par2par.dat.tpl
    %PEST_FILES_DIR%\pest_groups.txt
    %PEST_FILES_DIR%\pest_params.txt
    %PEST_FILES_DIR%\pest_webmod.bat.tpl
    %PEST_FILES_DIR%\tsproc.dat.tpl
    %PEST_FILES_DIR%\tsproc.in
    %PEST_FILES_DIR%\webmod.hydro.obs.ssf
    %PEST_FILES_DIR%\webmod.params.tpl
    %PEST_FILES_DIR%\webmod.statvar
    ) do (
    if NOT exist %%i (
        echo Did not find file %%i
        set errors=2
    )
)

REM check chemistry files
for %%i in ( %INPUT_DIR%\phreeq_lut 
    %INPUT_DIR%\phreeqc_web_lite.dat
    %INPUT_DIR%\webmod.chem.dat
    %INPUT_DIR%\webmod.pqi
    %PEST_FILES_DIR%\webmod.chem.obs.ssf
    ) do (
    if NOT exist %%i (
        echo Did not find chemistry file %%i
    )
)

REM quit if errors in hydro files
if %errors% == 2 (
    echo Stopping.
    goto :FINISH_UP
)
    
REM set port number
set PORT=4004 

REM setup output directory 
if NOT exist %OUTPUT_DIR% mkdir %OUTPUT_DIR%

REM setup working directory PROJECT_DIR=pest_run_dir
if exist %PROJECT_DIR% rmdir /s/q %PROJECT_DIR%
mkdir %PROJECT_DIR%

REM Sed to make tsproc.dat 
call :SED @PROJECT_DIR@ %PROJECT_DIR_SED%/ %PEST_FILES_DIR%\tsproc.dat.tpl > %PROJECT_DIR%\tsproc.dat

REM Sed to make pest_webmod.bat in PROJECT_DIR
call :SED @PROJECT_DIR@ %PROJECT_DIR_SED%/       %PEST_FILES_DIR%\pest_webmod.bat.tpl > %PROJECT_DIR%\pest_webmod.bat1
call :SED @PEST_BIN_DIR@ %PEST_BIN_DIR_SED%/     %PROJECT_DIR%\pest_webmod.bat1 >       %PROJECT_DIR%\pest_webmod.bat2
call :SED @TSPROC_BIN_DIR@ %TSPROC_BIN_DIR_SED%/ %PROJECT_DIR%\pest_webmod.bat2 >       %PROJECT_DIR%\pest_webmod.bat3
call :SED @DEL@ DEL                              %PROJECT_DIR%\pest_webmod.bat3 >       %PROJECT_DIR%\pest_webmod.bat
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
call :SED   ./input  .  %CONTROL_DIR%\webmod.control  >  %PROJECT_DIR%\webmod.control1
call :SED  ./output  .  %PROJECT_DIR%\webmod.control1 >  %PROJECT_DIR%\webmod.control2
call :SED ../../bin  .  %PROJECT_DIR%\webmod.control2 >  %PROJECT_DIR%\webmod.control
DEL %PROJECT_DIR%\webmod.control1  
DEL %PROJECT_DIR%\webmod.control2

REM change tsproc mode from "pest_prep" to "model_run"
call :SED  "CONTEXT pest_prep" "CONTEXT model_run" %PROJECT_DIR%\tsproc.dat > %PROJECT_DIR%\xxx
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
        ECHO 10.0  -3.0    0.3    0.03     -%NODES%  999  LAMFORGIVE >> %tmppest%
    ) else if "!COUNTER!" == "7" ( 
        ECHO 0.2   2.0   1.0e-3 >> %tmppest%
    ) else if "!COUNTER!" == "9" ( 
        ECHO 30   .005  4   4  .005   4 >> %tmppest%
    ) else if "!COUNTER!" == "13" ( 
        ECHO 32 5e-7 >> %tmppest%
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
%PEST_BIN_DIR%\beopest64.exe %PROJECT_DIR%\%PST% /H /L /p1 :%PORT% 

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
    copy webmod.hydro.out   ..\pest_results
    copy webmod.chem.out    ..\pest_results
    copy webmod.statvar     ..\pest_results
    copy webmod.pqi         ..\input\webmod.pqi.pest
    copy webmod.params      ..\input\webmod.params.pest
    copy webmod.hydro.out   ..\output\webmod.hydro.out.pest
    copy webmod.chem.out    ..\output\webmod.chem.out.pest
    copy webmod.statvar     ..\output\webmod.statvar.pest
    copy tsproc.dat         ..\pest_results
    cd ..\pest_results
    mkdir pest_files
    copy ..\pest_files\par2par.dat.tpl        .\pest_files
    copy ..\pest_files\pest_groups.txt        .\pest_files
    copy ..\pest_files\pest_params.txt        .\pest_files
    copy ..\pest_files\pest_webmod.bat.tpl    .\pest_files
    copy ..\pest_files\tsproc.dat.tpl         .\pest_files
    copy ..\pest_files\tsproc.in              .\pest_files
    copy ..\pest_files\webmod.chem.obs.ssf    .\pest_files
    copy ..\pest_files\webmod.hydro.obs.ssf   .\pest_files
    copy ..\pest_files\webmod.params.tpl      .\pest_files
    copy ..\pest_files\webmod.pqi.tpl         .\pest_files
    copy ..\pest_files\webmod.statvar         .\pest_files
REM Plot residuals and correlation
    %PEST_BIN_DIR%\pest_plot webmod.rei webmod__fit.pdf none
REM Plot parameters changes
    %PEST_BIN_DIR%\parm_plot webmod.pst webmod.sen webmod_par_calib.pdf none
REM Plot sensitivities
    %PEST_BIN_DIR%\sen_plot webmod.sen webmod_sensitivity.pdf
REM Plot contributions to phi by observation group
    %PEST_BIN_DIR%\pcon_plot webmod.rec webmod_phi.pdf none
REM Compute influence statistics
    %PEST_BIN_DIR%\infstat webmod webmod.infstat
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
echo usage: webpest_pc n  
echo where  n is number of workers to start
GOTO :EOF

:SED
REM -- Prepare the Command Processor --
SETLOCAL ENABLEEXTENSIONS
SETLOCAL DISABLEDELAYEDEXPANSION

::BatchSubstitude - parses a File line by line and replaces a substring"
::syntax: BatchSubstitude.bat OldStr NewStr File
::          OldStr [in] - string to be replaced
::          NewStr [in] - string to replace with
::          File   [in] - file to be parsed
:$changed 20100115
:$source http://www.dostips.com
if "%~1"=="" findstr "^::" "%~f0"&GOTO:EOF
for /f "tokens=1,* delims=]" %%A in ('"type %3|find /n /v """') do (
    set "line=%%B"
    if defined line (
        call set "line=echo.%%line:%~1=%~2%%"
        for /f "delims=" %%X in ('"echo."%%line%%""') do %%~X
    ) ELSE echo.
)
EXIT /b

:FINISH_UP
cd %PROJECT_DIR%\..


