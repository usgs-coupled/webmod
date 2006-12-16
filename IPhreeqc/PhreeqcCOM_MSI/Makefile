# $Id: Makefile 1452 2006-11-04 04:18:32Z charlton $

PHREEQC_COM_DIR  := ../PhreeqcCOM

PHREEQC_COM_DLL  := $(PHREEQC_COM_DIR)/Release/PhreeqcCOM.dll
PHREEQC_COM_TLB  := $(PHREEQC_COM_DIR)/Release/PhreeqcCOM.tlb
PHREEQC_COM_CHM  := $(PHREEQC_COM_DIR)/help/PhreeqcCOM.chm

CANDLE = candle
CANDLE_FLAGS=\
	-nologo \
	-dMergeDir=$(MERGE_DIR) \
	-dPhreeqcCOM_DLL=$(PHREEQC_COM_DLL) \
	-dPhreeqcCOM_TLB=$(PHREEQC_COM_TLB) \
	-dPhreeqcCOM_CHM=$(PHREEQC_COM_CHM)

LIGHT = light
LIGHT_FLAGS=\
	-nologo

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

PhreeqcCOM.msi : PhreeqcCOM.wixobj regsvr.wixobj $(PHREEQC_COM_DLL) $(PHREEQC_COM_TLB) $(PHREEQC_COM_CHM) $(BINARIES)
	$(LIGHT) $(LIGHT_FLAGS) PhreeqcCOM.wixobj regsvr.wixobj -out PhreeqcCOM.msi

PhreeqcCOM.wixobj : PhreeqcCOM.wxs Version.wxs
	$(CANDLE) $(CANDLE_FLAGS) -out $@ PhreeqcCOM.wxs

regsvr.wixobj : regsvr.wxs
	$(CANDLE) $(CANDLE_FLAGS) -out $@ regsvr.wxs

install :
	$(MSIEXEC) /i PhreeqcCOM.msi /qn+ INSTALLLOCATION=C:\DOCUME~1\charlton\PROGRA~1\PhreeqcCOM /log PhreeqcCOM.log
	
uninstall :
	$(MSIEXEC) /x PhreeqcCOM.msi /qn+

check : PhreeqcCOM.msi
	$(MSIVAL2) PhreeqcCOM.msi 'C:\Program Files\MsiVal2\logo.cub' -f

