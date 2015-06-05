REM @echo off
REM Sychronize latest WEBMOD executable, create the tsproc and PEST command files, and start beopest as master
REM Help printed if command entered with no arguments
IF %1.==. GOTO usage
IF %2.==. GOTO NoCase

REM Set variables
set nodes=%1
set pst=%2
set "WEB_DAT=..\WEBDATA"
set "WEB_TPL=..\WEBTPL"
set "WEB_SRC=.\webfiles"

set "PROJECT_DIR_PATH=C:\Programs\webmod-trunk\pest\andcrk\PC_Beotest\webmod"
set "PROJECT_DIR=C:\\Programs\\webmod-trunk\\pest\\andcrk\\PC_Beotest\\webmod"
set "PEST_BIN_DIR_PATH=C:\Programs\webmod-trunk\pest\PC_bin"
set "PEST_BIN_DIR=C:\\Programs\\webmod-trunk\\pest\\PC_bin"

REM dir
set PORT=4004 

REM setup WEB_SRC
if exist %WEB_SRC% rmdir /s/q %WEB_SRC%
mkdir %WEB_SRC%
REM sed is on NIXE so use it here. VB scripts are an alternative
sed "s#@PROJECT_DIR@#%PROJECT_DIR%\\#" %WEB_TPL%\andcrk.control.tpl > %WEB_SRC%\andcrk.control
sed "s#@PROJECT_DIR@#%PROJECT_DIR%\\#" %WEB_TPL%\andcrk_tsproc.tpl > %WEB_SRC%\andcrk_tsproc.dat
sed "s#@PROJECT_DIR@#%PROJECT_DIR%\\#g" %WEB_TPL%\pest_webmod.bat.tpl > %WEB_SRC%\pest_webmod.bata
sed "s#@PEST_BIN_DIR@#%PEST_BIN_DIR%\\#" %WEB_SRC%\pest_webmod.bata > %WEB_SRC%\pest_webmod.bat
del %WEB_SRC%\pest_webmod.bata
sed "s#@PROJECT_DIR@#%PROJECT_DIR%\\#g" %WEB_TPL%\tsproc.in.tpl > %WEB_SRC%\tsproc.in
sed "s#@PROJECT_DIR@#%PROJECT_DIR%\\#g" %WEB_TPL%\par2par_andcrk.tpl.tpl > %WEB_SRC%\par2par_andcrk.tpl
copy %WEB_DAT%\* %WEB_SRC%
copy %PEST_BIN_DIR_PATH%\webmod_1.0.exe %WEB_SRC%

REM setup working directory PROJECT_DIR_PATH
if exist %PROJECT_DIR_PATH% rmdir /s/q %PROJECT_DIR_PATH%
mkdir %PROJECT_DIR_PATH%
copy %WEB_SRC%\andcrk.statvar.full %PROJECT_DIR_PATH%\andcrk.statvar
XCOPY /I %WEB_SRC% %PROJECT_DIR_PATH%

REM run tsproc
cd %PROJECT_DIR_PATH%
%PEST_BIN_DIR_PATH%\tsproc.exe < %PROJECT_DIR_PATH%\tsproc.in
cd ..


sed -i "2s/CONTEXT pest_prep/CONTEXT model_run/" %PROJECT_DIR_PATH%\andcrk_tsproc.dat
sed -i "1s/.*/andcrk_tsproc.dat/" %PROJECT_DIR_PATH%\tsproc.in

REM Line 6 of the pst file RLAMBDA1 .....	
sed -i "6s/.*/10.0  -3.0    0.3    0.03     -10  999  LAMFORGIVE/" %PROJECT_DIR_PATH%\andcrk.pst
REM Line 7 of the pst file RELPARMAX FACPARMAX FACORIG.....	
sed -i "7s/.*/0.2   2.0   1.0e-3/" %PROJECT_DIR_PATH%\andcrk.pst
REM Line 9 of the pst file. NOPTMAX Max # of optimizations. Default is 30, set to 0 for single run with phi contributions, 1 for sensitivities, or a small number to test PEST loops.
sed -i "9s/.*/2   .005  4   4  .005   4/" %PROJECT_DIR_PATH%\andcrk.pst
REM Line 13. SVD block MAXSING EIGTHRESH. Replace MAXSING the maximum number of adjustable variables (number of singlular valuess at which truncation occurs)
sed -i "13s/.*/19 5e-7/" %PROJECT_DIR_PATH%\andcrk.pst
REM Line 14 EIGWRITE. 0 if not using SVD output file
sed -i "14s/.*/1/" %PROJECT_DIR_PATH%\andcrk.pst

REM remove sed files
rm sed*

REM Run parallel pest Master
cd %PROJECT_DIR_PATH%
start "Master" cmd /k call %PEST_BIN_DIR_PATH%\beopest64 %PROJECT_DIR_PATH%\%pst% /H :%PORT%
cd ..

sleep 5

REM start "Workers" cmd /k call .\slaves\runslave.bat %nodes% %pst%
REM goto :EOF

REM Run parallel pest workers
set MASTER=igskahhwwsdpark
set WEBDIR=%PROJECT_DIR_PATH%\
set WORKER_DIR=%PROJECT_DIR_PATH%\

REM make tmp dirs where slaves will be activated.
for /l %%X in (1, 1, %nodes%) do (
if exist %WORKER_DIR%\tmpest%%X rmdir /s/q %WORKER_DIR%\tmpest%%X 
mkdir %WORKER_DIR%\tmpest%%X
REM XCOPY /I %WEBDIR% %PROJECT_DIR_PATH%\tmpest%%X
cd %WORKER_DIR%\tmpest%%X
copy %WEBDIR%andcrk.pst .\
REM copy %WEBDIR%andcrk_tsproc.dat .\
REM copy %WEBDIR%andcrk.par2par.tpl .\
REM copy %WEBDIR%andcrk.statvar .\
copy %WEBDIR%params_andcrk.tpl .\
copy %WEBDIR%pqi_andcrk.tpl .\
copy %WEBDIR%pest_webmod.bat .\
copy %WEBDIR%andcrk.dat .\
copy %WEBDIR%andcrk.dat.chemdat .\
copy %WEBDIR%phreeq_lut .\
copy %WEBDIR%phreeqc_web_lite.dat .\
copy %WEBDIR%andcrk.ins .\
copy %WEBDIR%andcrk_tsproc.dat .\
START /B %PEST_BIN_DIR_PATH%\beopest64 %WEBDIR%%pst% /H %MASTER%:%PORT%
cd %WORKER_DIR%\
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


