REM Define ARCHIVE to be a directory for the WEBMOD Archive
REM Defined EXPORTED to be the directory of a WEBMOD trunk SVN export
REM EXPORTED is copied to ARCHIVE and unused files are deleted
REM
set ARCHIVE=D:\System\Users\rmwebb\Desktop\WEBMOD_1.0
set EXPORTED=D:\System\Users\rmwebb\Desktop\WebmodExport
REM
RMDIR /S /Q %EXPORTED%
svn export http://internalbrr.cr.usgs.gov/svn_GW/webmod/trunk %EXPORTED%
RMDIR /S /Q %ARCHIVE%
mkdir %ARCHIVE%
xcopy %EXPORTED% %ARCHIVE% /S /Y
cd %ARCHIVE%
DEL MakeArchive.bat
RMDIR bin /s /q
RMDIR lib /s /q
RMDIR projects /s /q
RMDIR tsproc-trunk /s /q
cd IPhreeqcMMS\IPhreeqc
del IPhreeqc*.* jenkins*.* Makefile.am vs*.* config* bootstrap All* INSTALL
RMDIR unit  /s /q
RMDIR tests  /s /q
RMDIR testcpp  /s /q
RMDIR test5  /s /q
RMDIR test2  /s /q
RMDIR test_ieee  /s /q
RMDIR test  /s /q
RMDIR R  /s /q
RMDIR phreeqc3-examples  /s /q
RMDIR packages /s /q
RMDIR memory_leak_f  /s /q
RMDIR examples  /s /q
RMDIR build  /s /q
RMDIR all /s /q
cd src\phreeqcpp
del Makefile*
cd %ARCHIVE%
mkdir               SourceCodeAndExecutables
move CMakeLists.txt SourceCodeAndExecutables\.
move IPhreeqcMMS    SourceCodeAndExecutables\.
move mmf_c          SourceCodeAndExecutables\.
move Tests          SourceCodeAndExecutables\.
move webmod         SourceCodeAndExecutables\.
cd %ARCHIVE%\..
PAUSE
