# Microsoft Developer Studio Project File - Name="IPhreeqc" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=IPhreeqc - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "IPhreeqc.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "IPhreeqc.mak" CFG="IPhreeqc - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "IPhreeqc - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "IPhreeqc - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
F90=df.exe
RSC=rc.exe

!IF  "$(CFG)" == "IPhreeqc - Win32 Release"

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
# ADD F90 /compile_only /nologo /warn:nofileopt
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "include" /D "NDEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "SWIG_SHARED_OBJ" /YX /FD /c
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"lib\IPhreeqc.lib"

!ELSEIF  "$(CFG)" == "IPhreeqc - Win32 Debug"

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
# ADD F90 /check:bounds /compile_only /debug:full /nologo /traceback /warn:argument_checking /warn:nofileopt
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /W3 /Gm /GX /ZI /Od /I "include" /D "_DEBUG" /D "WIN32" /D "_MBCS" /D "_LIB" /D "SWIG_SHARED_OBJ" /YX /FD /GZ /c
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"lib\IPhreeqcd.lib"

!ENDIF 

# Begin Target

# Name "IPhreeqc - Win32 Release"
# Name "IPhreeqc - Win32 Debug"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat;f90;for;f;fpp"
# Begin Group "phreeqc"

# PROP Default_Filter ""
# Begin Source File

SOURCE=.\src\phreeqc\advection.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\basic.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\basicsubs.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\cl1.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\cvdense.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\cvode.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\dense.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\integrate.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\inverse.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\isotopes.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\kinetics.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\main.c
# PROP Exclude_From_Build 1
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\mainsubs.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\model.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\nvector.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\nvector_serial.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\p2clib.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\parse.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\phqalloc.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\prep.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\print.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\read.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\readtr.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\smalldense.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\spread.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\step.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\structures.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\sundialsmath.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\tidy.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\transport.c
# End Source File
# Begin Source File

SOURCE=.\src\phreeqc\utilities.c
# End Source File
# End Group
# Begin Source File

SOURCE=.\src\fwrap.c
# End Source File
# Begin Source File

SOURCE=.\src\global.c
# End Source File
# Begin Source File

SOURCE=.\src\IPhreeqc.cxx
# End Source File
# Begin Source File

SOURCE=.\src\IPhreeqcF.f
# ADD F90 /fpp
# End Source File
# Begin Source File

SOURCE=.\src\Output.cxx
# End Source File
# Begin Source File

SOURCE=.\src\Overrides.cxx
# End Source File
# Begin Source File

SOURCE=.\src\Phreeqc.cxx
# End Source File
# Begin Source File

SOURCE=.\src\PhreeqcParser.cxx
# End Source File
# Begin Source File

SOURCE=.\src\SelectedOutput.cxx
# End Source File
# Begin Source File

SOURCE=.\src\Var.c
# End Source File
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl;fi;fd"
# Begin Source File

SOURCE=.\src\CVar.hxx
# End Source File
# Begin Source File

SOURCE=.\src\ErrorReporter.hxx
# End Source File
# Begin Source File

SOURCE=.\include\IPhreeqc.f90.inc
# End Source File
# Begin Source File

SOURCE=.\include\IPhreeqc.h
# End Source File
# Begin Source File

SOURCE=.\src\Output.hxx
# End Source File
# Begin Source File

SOURCE=.\src\Phreeqc.hxx
# End Source File
# Begin Source File

SOURCE=.\src\phreeqcns.hxx
# End Source File
# Begin Source File

SOURCE=.\src\PhreeqcParser.hxx
# End Source File
# Begin Source File

SOURCE=.\src\SelectedOutput.hxx
# End Source File
# Begin Source File

SOURCE=.\src\Var.h
# End Source File
# End Group
# End Target
# End Project
