set ARCHIVE=C:\Programs\webmod-trunk\webmod1.0
set EXPORTED=C:\Programs\webmod-trunk\webmodexport
DEL C:\Programs\webmod-trunk\webmod1.0 /Q /F
mkdir C:\Programs\webmod-trunk\webmod1.0
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
cd %ARCHIVE%\..
