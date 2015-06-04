REM @echo off
REM Sychronize latest WEBMOD executable, create the tsproc and PEST command files, and start beopest as master
REM Help printed if command entered with no arguments
IF %1.==. GOTO usage
REM 
set "WEB_DAT=..\WEBDATA"
set "WEB_TPL=..\WEBTPL"
set "WEB_SRC=.\webfiles"

set "PROJECT_DIR_PATH=C:\Programs\webmod-trunk\pest\andcrk\PC_Beotest\webmod"
set "PROJECT_DIR=C:\\Programs\\webmod-trunk\\pest\\andcrk\\PC_Beotest\\webmod"
set "PEST_BIN_DIR_PATH=C:\Programs\webmod-trunk\pest\PC_bin"
set "PEST_BIN_DIR=C:\\Programs\\webmod-trunk\\pest\\PC_bin"
REM dir
set PORT=4004 
set pst=%1
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
rename %WEB_SRC%\andcrk.statvar.full andcrk.statvar
if exist %PROJECT_DIR_PATH% rmdir /s/q %PROJECT_DIR_PATH%
copy %PEST_BIN_DIR_PATH%\webmod_1.0.exe %WEB_SRC%
XCOPY /I %WEB_SRC% %PROJECT_DIR_PATH%
copy %PROJECT_DIR_PATH%\andcrk.statvar .\
%PEST_BIN_DIR_PATH%\tsproc.exe < %PROJECT_DIR_PATH%\tsproc.in
copy andcrk.ins %PROJECT_DIR_PATH%

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

REM cd x:\webmod
REM x:
REM pwd
REM pause
REM ..\bin\beopest64 %PROJECT_DIR_PATH%\%pst% /H :%PORT%
%PEST_BIN_DIR_PATH%\beopest64 %PROJECT_DIR_PATH%\%pst% /H :%PORT%
GOTO :EOF

:usage
echo.
echo usage: runmaster case.pst
echo where  case.pst is the name of the PEST control file.
:EOF


