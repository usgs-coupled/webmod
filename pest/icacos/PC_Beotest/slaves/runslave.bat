@echo off
REM Help printed if command entered with no arguments
IF %1.==. GOTO usage
IF %2.==. GOTO NoCase
REM 
REM set WEBDIR=\\igskahhwuwhopkn\pest\webmod\
REM set PESTBINDIR=\\igskahhwuwhopkn\pest\bin\
REM set MASTER=igskahhwwsnixe

set WEBDIR=C:\Programs\webmod-trunk\pest\andcrk\PC_Beotest\webmod\
set PESTBINDIR=C:\Programs\webmod-trunk\pest\PC_bin\
set MASTER=igskahhwwsdpark

set PORT=4004
set nodes=%1
set pst=%2
REM make tmp dirs where slaves will be activated.
for /l %%X in (1, 1, %nodes%) do (
if exist tmpest%%X rmdir /s/q tmpest%%X 
mkdir tmpest%%X
REM XCOPY /I %WEBDIR% tmpest%%X
cd tmpest%%X
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
START /B %PESTBINDIR%beopest64 %WEBDIR%%pst% /H %MASTER%:%PORT%
cd ..
)
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
GOTO :EOF
