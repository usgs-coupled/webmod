@echo off
REM Sychronize latest WEBMOD executable, create the tsproc and PEST command files, and start beopest as master
REM Help printed if command entered with no arguments
IF %1.==. GOTO usage
IF %2.==. GOTO NoCase

set HOME_BATCH=%cd%
set HOME=%HOME_BATCH:\=/%

REM Set variables
set nodes=%1
set pst=%2
set BASIN=andcrk
set "WEB_DAT=%HOME_BATCH%\..\WEBDATA"
set "WEB_TPL=%HOME_BATCH%\..\WEBTPL"
set "WEB_SRC=.\webfiles"
set "PROJECT_DIR_PATH=%HOME_BATCH%\webmod"
set "PROJECT_DIR=%HOME%/webmod"
set "PEST_BIN_DIR_PATH=%HOME_BATCH%\..\..\PC_bin"
set "PEST_BIN_DIR=%HOME%/../../PC_bin"
set "TSPROC_BIN_DIR=%HOME%/../../PC_bin"

REM set port number
set PORT=4004 

REM setup working directory PROJECT_DIR_PATH=webmod
if exist %PROJECT_DIR_PATH% rmdir /s/q %PROJECT_DIR_PATH%
mkdir %PROJECT_DIR_PATH%

REM Sed to make basin_tsproc.dat 
sed "s#@PROJECT_DIR@#%PROJECT_DIR%/#" %WEB_TPL%\%BASIN%_tsproc.tpl > %PROJECT_DIR%\%BASIN%_tsproc.dat
sed -i "s#@PROJECT_DIR@#%PROJECT_DIR%/#" %PROJECT_DIR%\%BASIN%_tsproc.dat

REM Sed to make pest_webmod.bat in PROJECT_DIR
sed "s#@PROJECT_DIR@#%PROJECT_DIR%/#g" %WEB_TPL%\pest_webmod.bat.tpl > %PROJECT_DIR_PATH%\pest_webmod.bat
sed -i "s#@PEST_BIN_DIR@#%PEST_BIN_DIR%/#g" %PROJECT_DIR_PATH%\pest_webmod.bat
sed -i "s#@TSPROC_BIN_DIR@#%TSPROC_BIN_DIR%/#" %PROJECT_DIR_PATH%\pest_webmod.bat
sed -i "s#@DEL@#DEL#" %PROJECT_DIR_PATH%\pest_webmod.bat
REM remove sed files
rm sed*

REM Copy files to tsproc directory
copy %WEB_DAT%\*.ssf            	%PROJECT_DIR_PATH%
copy %WEB_DAT%\pest_groups*.txt 	%PROJECT_DIR_PATH%
copy %WEB_DAT%\pest_params*.txt 	%PROJECT_DIR_PATH%
copy %WEB_DAT%\tsproc.in        	%PROJECT_DIR_PATH%
copy %WEB_DAT%\%BASIN%.statvar.full	%PROJECT_DIR_PATH%\%BASIN%.statvar
copy %WEB_DAT%\par2par_%BASIN%.tpl	%PROJECT_DIR_PATH%

REM run tsproc
cd %PROJECT_DIR_PATH%
%PEST_BIN_DIR_PATH%\tsproc.exe < .\tsproc.in

REM move output from tsproc to working directory
cp %TSPROC_BIN_DIR%\webmod_1.0.exe  %PROJECT_DIR_PATH%
cp %WEB_DAT%\%BASIN%.control         %PROJECT_DIR_PATH%

REM change tsproc mode from "pest_prep" to "model_run"
sed -i "2s/CONTEXT pest_prep/CONTEXT model_run/" %PROJECT_DIR_PATH%\%BASIN%_tsproc.dat
REM Line 6 of the pst file RLAMBDA1 .....	
sed -i "6s/.*/10.0  -3.0    0.3    0.03     -10  999  LAMFORGIVE/" %PROJECT_DIR_PATH%\%BASIN%.pst
REM Line 7 of the pst file RELPARMAX FACPARMAX FACORIG.....	
sed -i "7s/.*/0.2   2.0   1.0e-3/" %PROJECT_DIR_PATH%\%BASIN%.pst
REM Line 9 of the pst file. NOPTMAX Max # of optimizations. Default is 30, set to 0 for single run with phi contributions, 1 for sensitivities, or a small number to test PEST loops.
sed -i "9s/.*/30   .005  4   4  .005   4/" %PROJECT_DIR_PATH%\%BASIN%.pst
REM Line 13. SVD block MAXSING EIGTHRESH. Replace MAXSING the maximum number of adjustable variables (number of singlular valuess at which truncation occurs)
sed -i "13s/.*/19 5e-7/" %PROJECT_DIR_PATH%\%BASIN%.pst
REM Line 14 EIGWRITE. 0 if not using SVD output file
sed -i "14s/.*/1/" %PROJECT_DIR_PATH%\%BASIN%.pst
REM remove sed files
rm sed*
cd %PROJECT_DIR_PATH%

REM Run parallel pest Master
cd %PROJECT_DIR_PATH%
start "Master" cmd /k "call %PEST_BIN_DIR_PATH%\beopest64 %PROJECT_DIR_PATH%\%pst% /H :%PORT% & cd .."

REM wait 2 seconds for master to come up before starting workers
sleep 2

REM make tmp dirs and run parallel pest workers
set MASTER=%COMPUTERNAME%
for /l %%X in (1, 1, %nodes%) do (
	if exist %PROJECT_DIR_PATH%\tmpest%%X rmdir /s/q %PROJECT_DIR_PATH%\tmpest%%X 
	mkdir %PROJECT_DIR_PATH%\tmpest%%X
	cd %PROJECT_DIR_PATH%\tmpest%%X
	copy %WEB_DAT%\params_%BASIN%.tpl .\
	copy %WEB_DAT%\pqi_%BASIN%.tpl .\
	copy %WEB_DAT%\%BASIN%.dat .\
	copy %WEB_DAT%\%BASIN%.dat.chemdat .\
	copy %WEB_DAT%\phreeq_lut .\
	copy %WEB_DAT%\phreeqc_web_lite.dat .\
	copy %PROJECT_DIR_PATH%\pest_webmod.bat .\
	copy %PROJECT_DIR_PATH%\%BASIN%.ins .\
	copy %PROJECT_DIR_PATH%\%BASIN%_tsproc.dat .\
	START /B %PEST_BIN_DIR_PATH%\beopest64 %PROJECT_DIR_PATH%\%pst% /H %MASTER%:%PORT%
	cd %PROJECT_DIR_PATH%
)
cd %PROJECT_DIR_PATH%\..
GOTO :EOF

:usage
echo.
echo usage: runslave n case.pst
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


