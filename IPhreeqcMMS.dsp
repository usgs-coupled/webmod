# Microsoft Developer Studio Project File - Name="IPhreeqcMMS" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=IPhreeqcMMS - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "IPhreeqcMMS.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "IPhreeqcMMS.mak" CFG="IPhreeqcMMS - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "IPhreeqcMMS - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "IPhreeqcMMS - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "IPhreeqcMMS - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE F90 /compile_only /nologo /warn:nofileopt
# ADD F90 /compile_only /include:"IPhreeqc/include" /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "$(DEV_GMP_INC)" /I "include" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "SWIG_SHARED_OBJ" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"lib\IPhreeqcMMS.lib"

!ELSEIF  "$(CFG)" == "IPhreeqcMMS - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD F90 /browser /check:bounds /compile_only /debug:full /include:"IPhreeqc/include" /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "$(DEV_GMP_INC)" /I "include" /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "SWIG_SHARED_OBJ" /FR /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"lib\IPhreeqcMMSd.lib"

!ENDIF 

# Begin Target

# Name "IPhreeqcMMS - Win32 Release"
# Name "IPhreeqcMMS - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Group "IPhreeqc"

# PROP Default_Filter ""
# Begin Group "phreeqcpp"

# PROP Default_Filter ""
# Begin Group "phreeqc"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\advection.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\basic.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\basicsubs.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\cl1.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\cl1mp.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\cvdense.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\cvode.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\dense.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\dw.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\input.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\integrate.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\inverse.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\isotopes.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\kinetics.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\main.c
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\mainsubs.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\model.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\nvector.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\nvector_serial.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\output.c
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\p2clib.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\parse.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\phqalloc.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\phreeqc_files.c
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\pitzer.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\pitzer_structures.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\prep.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\print.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\read.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\readtr.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\smalldense.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\spread.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\step.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\structures.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\sundialsmath.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\tally.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\tidy.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\transport.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcpp\phreeqc\utilities.c
# End Source File
# End Group
# End Group
# Begin Source File

SOURCE=.\IPhreeqc\src\fwrap.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\global.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\IPhreeqc.cxx
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\module_files.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\module_output.c
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\Phreeqc.cxx
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\SelectedOutput.cxx
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\TallyF.F
# ADD F90 /fpp
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\Var.c
# End Source File
# End Group
# Begin Source File

SOURCE=.\src\phr_cmix.c
# End Source File
# Begin Source File

SOURCE=.\src\phr_mix.F
DEP_F90_PHR_M=\
	".\IPhreeqc\include\IPhreeqc.f.inc"\
	
# ADD F90 /fpp
# End Source File
# Begin Source File

SOURCE=.\src\phr_multicopy.f
DEP_F90_PHR_MU=\
	".\IPhreeqc\include\IPhreeqc.f.inc"\
	
# ADD F90 /include:"../IPhreeqc/include/"
# End Source File
# Begin Source File

SOURCE=.\src\phr_precip.f
DEP_F90_PHR_P=\
	".\IPhreeqc\include\IPhreeqc.f.inc"\
	
# ADD F90 /include:"../IPhreeqc/include/"
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\IPhreeqc\src\CVar.hxx
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\ErrorReporter.hxx
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\include\IPhreeqc.h
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\Output.hxx
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\Phreeqc.hxx
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\phreeqcns.hxx
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\PhreeqcParser.hxx
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\src\SelectedOutput.hxx
# End Source File
# Begin Source File

SOURCE=.\IPhreeqc\include\Var.h
# End Source File
# End Group
# End Target
# End Project
