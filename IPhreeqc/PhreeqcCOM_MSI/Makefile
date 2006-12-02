# $Id: Makefile 1452 2006-11-04 04:18:32Z charlton $

DEVENV := /cygdrive/c/Program\ Files/Microsoft\ Visual\ Studio\ 8/Common7/IDE/devenv

SLN := ../WPhast.sln

PHREEQC_COM_DIR := ../PhreeqcCOM/Release
PHREEQC_COM := $(PHREEQC_COM_DIR)/PhreeqcCOM.dll

CANDLE = candle
CANDLE_FLAGS=\
	-nologo \
	-dMergeDir=$(MERGE_DIR) \
	-dPhreeqcCOMDir=$(PHREEQC_COM_DIR)

LIGHT = light
LIGHT_FLAGS=\
	-nologo

CSCRIPT = cscript

MERGE_DIR := "C:\Program Files\Common Files\Merge Modules"

MSIVAL2 := /cygdrive/c/Program\ Files/MsiVal2/msival2

MSIEXEC := msiexec

BINARIES=\
	Binary/Banner.bmp\
	Binary/Banner.bmp\
	Binary/Complete.ico\
	Binary/Custom.ico\
	Binary/Dialog.bmp\
	Binary/Exclam.ico\
	Binary/Info.ico\
	Binary/License.rtf\
	Binary/New.ico\
	Binary/Remove.ico\
	Binary/Repair.ico\
	Binary/Typical.ico\
	Binary/Up.ico
	
all : PhreeqcCOM.msi

PhreeqcCOM.msi : PhreeqcCOM.wixobj regsvr.wixobj $(PHREEQC_COM) $(BINARIES)
	$(LIGHT) $(LIGHT_FLAGS) PhreeqcCOM.wixobj regsvr.wixobj -out PhreeqcCOM.msi

PhreeqcCOM.wixobj : PhreeqcCOM.wxs Version.wxs
	$(CANDLE) $(CANDLE_FLAGS) -out $@ PhreeqcCOM.wxs

regsvr.wixobj : regsvr.wxs
	$(CANDLE) $(CANDLE_FLAGS) -out $@ regsvr.wxs

##$(WPHAST_EXE) : $(SLN)
##	@rm -rf Release.log
##	$(DEVENV) /out Release.log /build Release $(SLN)
##	@cat Release.log
	
install :
	$(MSIEXEC) /i PhreeqcCOM.msi /qn+ INSTALLLOCATION=C:\DOCUME~1\charlton\PROGRA~1\PhreeqcCOM /log PhreeqcCOM.log
	
uninstall :
	$(MSIEXEC) /x PhreeqcCOM.msi /qn+

check : PhreeqcCOM.msi
	$(MSIVAL2) PhreeqcCOM.msi 'C:\Program Files\MsiVal2\logo.cub' -f

