# Compiler flags...
CPP_COMPILER = g++
C_COMPILER = gcc

# Include paths...
Debug_Include_Path=-I"src" -I"IPhreeqc/src" -I"IPhreeqc/src/phreeqcpp" 
Debug_Include_Path=-I"src" -I"IPhreeqc/src" -I"IPhreeqc/src/phreeqcpp" 
Release_Include_Path=-I"src" -I"IPhreeqc/src" -I"IPhreeqc/src/phreeqcpp" 
Release_Include_Path=-I"src" -I"IPhreeqc/src" -I"IPhreeqc/src/phreeqcpp" 
Template_Include_Path=-I"src" -I"IPhreeqc/src" -I"IPhreeqc/src/phreeqcpp" 
Template_Include_Path=-I"src" -I"IPhreeqc/src" -I"IPhreeqc/src/phreeqcpp" 

# Library paths...
Debug_Library_Path=
Debug_Library_Path=
Release_Library_Path=
Release_Library_Path=
Template_Library_Path=
Template_Library_Path=

# Additional libraries...
Debug_Libraries=
Debug_Libraries=
Release_Libraries=
Release_Libraries=
Template_Libraries=
Template_Libraries=

# Preprocessor definitions...
Debug_Preprocessor_Definitions=-D _DEBUG -D GCC_BUILD -D _LIB -D SWIG_SHARED_OBJ -D _CRT_SECURE_NO_DEPRECATE 
Debug_Preprocessor_Definitions=-D _DEBUG -D GCC_BUILD -D _LIB -D SWIG_SHARED_OBJ -D _CRT_SECURE_NO_DEPRECATE 
Release_Preprocessor_Definitions=-D NDEBUG -D GCC_BUILD -D _LIB -D SWIG_SHARED_OBJ -D _CRT_SECURE_NO_DEPRECATE 
Release_Preprocessor_Definitions=-D NDEBUG -D GCC_BUILD -D _LIB -D SWIG_SHARED_OBJ -D _CRT_SECURE_NO_DEPRECATE 
Template_Preprocessor_Definitions=-D GCC_BUILD 
Template_Preprocessor_Definitions=-D GCC_BUILD 

# Implictly linked object files...
Debug_Implicitly_Linked_Objects=
Debug_Implicitly_Linked_Objects=
Release_Implicitly_Linked_Objects=
Release_Implicitly_Linked_Objects=
Template_Implicitly_Linked_Objects=
Template_Implicitly_Linked_Objects=

# Compiler flags...
Debug_Compiler_Flags=-O0 
Debug_Compiler_Flags=-O0 
Release_Compiler_Flags=-O2 
Release_Compiler_Flags=-O2 
Template_Compiler_Flags=-O2 
Template_Compiler_Flags=-O2 

# Builds all configurations for this project...
.PHONY: build_all_configurations
build_all_configurations: Debug Debug Release Release Template Template 

# Builds the Debug configuration...
.PHONY: Debug
Debug: create_folders gccDebug/src/cdecl.o gccDebug/src/fortran.o gccDebug/src/IPhreeqcMMS.o gccDebug/src/IPhreeqcMMSLib.o gccDebug/src/stdcall.o 
	g++ gccDebug/src/cdecl.o gccDebug/src/fortran.o gccDebug/src/IPhreeqcMMS.o gccDebug/src/IPhreeqcMMSLib.o gccDebug/src/stdcall.o  $(Debug_Library_Path) $(Debug_Libraries) -Wl,-rpath,./ -o ../gccDebug/IPhreeqcMMS.exe

# Compiles file src/cdecl.cpp for the Debug configuration...
-include gccDebug/src/cdecl.d
gccDebug/src/cdecl.o: src/cdecl.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/cdecl.cpp $(Debug_Include_Path) -o gccDebug/src/cdecl.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/cdecl.cpp $(Debug_Include_Path) > gccDebug/src/cdecl.d

# Compiles file src/fortran.cpp for the Debug configuration...
-include gccDebug/src/fortran.d
gccDebug/src/fortran.o: src/fortran.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/fortran.cpp $(Debug_Include_Path) -o gccDebug/src/fortran.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/fortran.cpp $(Debug_Include_Path) > gccDebug/src/fortran.d

# Compiles file src/IPhreeqcMMS.cpp for the Debug configuration...
-include gccDebug/src/IPhreeqcMMS.d
gccDebug/src/IPhreeqcMMS.o: src/IPhreeqcMMS.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/IPhreeqcMMS.cpp $(Debug_Include_Path) -o gccDebug/src/IPhreeqcMMS.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/IPhreeqcMMS.cpp $(Debug_Include_Path) > gccDebug/src/IPhreeqcMMS.d

# Compiles file src/IPhreeqcMMSLib.cpp for the Debug configuration...
-include gccDebug/src/IPhreeqcMMSLib.d
gccDebug/src/IPhreeqcMMSLib.o: src/IPhreeqcMMSLib.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/IPhreeqcMMSLib.cpp $(Debug_Include_Path) -o gccDebug/src/IPhreeqcMMSLib.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/IPhreeqcMMSLib.cpp $(Debug_Include_Path) > gccDebug/src/IPhreeqcMMSLib.d

# Compiles file src/stdcall.cpp for the Debug configuration...
-include gccDebug/src/stdcall.d
gccDebug/src/stdcall.o: src/stdcall.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/stdcall.cpp $(Debug_Include_Path) -o gccDebug/src/stdcall.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/stdcall.cpp $(Debug_Include_Path) > gccDebug/src/stdcall.d

# Builds the Debug configuration...
.PHONY: Debug
Debug: create_folders x64/gccDebug/src/cdecl.o x64/gccDebug/src/fortran.o x64/gccDebug/src/IPhreeqcMMS.o x64/gccDebug/src/IPhreeqcMMSLib.o x64/gccDebug/src/stdcall.o 
	g++ x64/gccDebug/src/cdecl.o x64/gccDebug/src/fortran.o x64/gccDebug/src/IPhreeqcMMS.o x64/gccDebug/src/IPhreeqcMMSLib.o x64/gccDebug/src/stdcall.o  $(Debug_Library_Path) $(Debug_Libraries) -Wl,-rpath,./ -o ../x64/gccDebug/IPhreeqcMMS.exe

# Compiles file src/cdecl.cpp for the Debug configuration...
-include x64/gccDebug/src/cdecl.d
x64/gccDebug/src/cdecl.o: src/cdecl.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/cdecl.cpp $(Debug_Include_Path) -o x64/gccDebug/src/cdecl.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/cdecl.cpp $(Debug_Include_Path) > x64/gccDebug/src/cdecl.d

# Compiles file src/fortran.cpp for the Debug configuration...
-include x64/gccDebug/src/fortran.d
x64/gccDebug/src/fortran.o: src/fortran.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/fortran.cpp $(Debug_Include_Path) -o x64/gccDebug/src/fortran.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/fortran.cpp $(Debug_Include_Path) > x64/gccDebug/src/fortran.d

# Compiles file src/IPhreeqcMMS.cpp for the Debug configuration...
-include x64/gccDebug/src/IPhreeqcMMS.d
x64/gccDebug/src/IPhreeqcMMS.o: src/IPhreeqcMMS.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/IPhreeqcMMS.cpp $(Debug_Include_Path) -o x64/gccDebug/src/IPhreeqcMMS.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/IPhreeqcMMS.cpp $(Debug_Include_Path) > x64/gccDebug/src/IPhreeqcMMS.d

# Compiles file src/IPhreeqcMMSLib.cpp for the Debug configuration...
-include x64/gccDebug/src/IPhreeqcMMSLib.d
x64/gccDebug/src/IPhreeqcMMSLib.o: src/IPhreeqcMMSLib.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/IPhreeqcMMSLib.cpp $(Debug_Include_Path) -o x64/gccDebug/src/IPhreeqcMMSLib.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/IPhreeqcMMSLib.cpp $(Debug_Include_Path) > x64/gccDebug/src/IPhreeqcMMSLib.d

# Compiles file src/stdcall.cpp for the Debug configuration...
-include x64/gccDebug/src/stdcall.d
x64/gccDebug/src/stdcall.o: src/stdcall.cpp
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/stdcall.cpp $(Debug_Include_Path) -o x64/gccDebug/src/stdcall.o
	$(CPP_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -MM src/stdcall.cpp $(Debug_Include_Path) > x64/gccDebug/src/stdcall.d

# Builds the Release configuration...
.PHONY: Release
Release: create_folders gccRelease/src/cdecl.o gccRelease/src/fortran.o gccRelease/src/IPhreeqcMMS.o gccRelease/src/IPhreeqcMMSLib.o gccRelease/src/stdcall.o 
	g++ gccRelease/src/cdecl.o gccRelease/src/fortran.o gccRelease/src/IPhreeqcMMS.o gccRelease/src/IPhreeqcMMSLib.o gccRelease/src/stdcall.o  $(Release_Library_Path) $(Release_Libraries) -Wl,-rpath,./ -o ../gccRelease/IPhreeqcMMS.exe

# Compiles file src/cdecl.cpp for the Release configuration...
-include gccRelease/src/cdecl.d
gccRelease/src/cdecl.o: src/cdecl.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/cdecl.cpp $(Release_Include_Path) -o gccRelease/src/cdecl.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/cdecl.cpp $(Release_Include_Path) > gccRelease/src/cdecl.d

# Compiles file src/fortran.cpp for the Release configuration...
-include gccRelease/src/fortran.d
gccRelease/src/fortran.o: src/fortran.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/fortran.cpp $(Release_Include_Path) -o gccRelease/src/fortran.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/fortran.cpp $(Release_Include_Path) > gccRelease/src/fortran.d

# Compiles file src/IPhreeqcMMS.cpp for the Release configuration...
-include gccRelease/src/IPhreeqcMMS.d
gccRelease/src/IPhreeqcMMS.o: src/IPhreeqcMMS.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/IPhreeqcMMS.cpp $(Release_Include_Path) -o gccRelease/src/IPhreeqcMMS.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/IPhreeqcMMS.cpp $(Release_Include_Path) > gccRelease/src/IPhreeqcMMS.d

# Compiles file src/IPhreeqcMMSLib.cpp for the Release configuration...
-include gccRelease/src/IPhreeqcMMSLib.d
gccRelease/src/IPhreeqcMMSLib.o: src/IPhreeqcMMSLib.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/IPhreeqcMMSLib.cpp $(Release_Include_Path) -o gccRelease/src/IPhreeqcMMSLib.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/IPhreeqcMMSLib.cpp $(Release_Include_Path) > gccRelease/src/IPhreeqcMMSLib.d

# Compiles file src/stdcall.cpp for the Release configuration...
-include gccRelease/src/stdcall.d
gccRelease/src/stdcall.o: src/stdcall.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/stdcall.cpp $(Release_Include_Path) -o gccRelease/src/stdcall.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/stdcall.cpp $(Release_Include_Path) > gccRelease/src/stdcall.d

# Builds the Release configuration...
.PHONY: Release
Release: create_folders x64/gccRelease/src/cdecl.o x64/gccRelease/src/fortran.o x64/gccRelease/src/IPhreeqcMMS.o x64/gccRelease/src/IPhreeqcMMSLib.o x64/gccRelease/src/stdcall.o 
	g++ x64/gccRelease/src/cdecl.o x64/gccRelease/src/fortran.o x64/gccRelease/src/IPhreeqcMMS.o x64/gccRelease/src/IPhreeqcMMSLib.o x64/gccRelease/src/stdcall.o  $(Release_Library_Path) $(Release_Libraries) -Wl,-rpath,./ -o ../x64/gccRelease/IPhreeqcMMS.exe

# Compiles file src/cdecl.cpp for the Release configuration...
-include x64/gccRelease/src/cdecl.d
x64/gccRelease/src/cdecl.o: src/cdecl.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/cdecl.cpp $(Release_Include_Path) -o x64/gccRelease/src/cdecl.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/cdecl.cpp $(Release_Include_Path) > x64/gccRelease/src/cdecl.d

# Compiles file src/fortran.cpp for the Release configuration...
-include x64/gccRelease/src/fortran.d
x64/gccRelease/src/fortran.o: src/fortran.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/fortran.cpp $(Release_Include_Path) -o x64/gccRelease/src/fortran.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/fortran.cpp $(Release_Include_Path) > x64/gccRelease/src/fortran.d

# Compiles file src/IPhreeqcMMS.cpp for the Release configuration...
-include x64/gccRelease/src/IPhreeqcMMS.d
x64/gccRelease/src/IPhreeqcMMS.o: src/IPhreeqcMMS.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/IPhreeqcMMS.cpp $(Release_Include_Path) -o x64/gccRelease/src/IPhreeqcMMS.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/IPhreeqcMMS.cpp $(Release_Include_Path) > x64/gccRelease/src/IPhreeqcMMS.d

# Compiles file src/IPhreeqcMMSLib.cpp for the Release configuration...
-include x64/gccRelease/src/IPhreeqcMMSLib.d
x64/gccRelease/src/IPhreeqcMMSLib.o: src/IPhreeqcMMSLib.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/IPhreeqcMMSLib.cpp $(Release_Include_Path) -o x64/gccRelease/src/IPhreeqcMMSLib.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/IPhreeqcMMSLib.cpp $(Release_Include_Path) > x64/gccRelease/src/IPhreeqcMMSLib.d

# Compiles file src/stdcall.cpp for the Release configuration...
-include x64/gccRelease/src/stdcall.d
x64/gccRelease/src/stdcall.o: src/stdcall.cpp
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/stdcall.cpp $(Release_Include_Path) -o x64/gccRelease/src/stdcall.o
	$(CPP_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -MM src/stdcall.cpp $(Release_Include_Path) > x64/gccRelease/src/stdcall.d

# Builds the Template configuration...
.PHONY: Template
Template: create_folders gccTemplate/src/cdecl.o gccTemplate/src/fortran.o gccTemplate/src/IPhreeqcMMS.o gccTemplate/src/IPhreeqcMMSLib.o gccTemplate/src/stdcall.o 
	g++ gccTemplate/src/cdecl.o gccTemplate/src/fortran.o gccTemplate/src/IPhreeqcMMS.o gccTemplate/src/IPhreeqcMMSLib.o gccTemplate/src/stdcall.o  $(Template_Library_Path) $(Template_Libraries) -Wl,-rpath,./ -o ../gccTemplate/IPhreeqcMMS.exe

# Compiles file src/cdecl.cpp for the Template configuration...
-include gccTemplate/src/cdecl.d
gccTemplate/src/cdecl.o: src/cdecl.cpp
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -c src/cdecl.cpp $(Template_Include_Path) -o gccTemplate/src/cdecl.o
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -MM src/cdecl.cpp $(Template_Include_Path) > gccTemplate/src/cdecl.d

# Compiles file src/fortran.cpp for the Template configuration...
-include gccTemplate/src/fortran.d
gccTemplate/src/fortran.o: src/fortran.cpp
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -c src/fortran.cpp $(Template_Include_Path) -o gccTemplate/src/fortran.o
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -MM src/fortran.cpp $(Template_Include_Path) > gccTemplate/src/fortran.d

# Compiles file src/IPhreeqcMMS.cpp for the Template configuration...
-include gccTemplate/src/IPhreeqcMMS.d
gccTemplate/src/IPhreeqcMMS.o: src/IPhreeqcMMS.cpp
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -c src/IPhreeqcMMS.cpp $(Template_Include_Path) -o gccTemplate/src/IPhreeqcMMS.o
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -MM src/IPhreeqcMMS.cpp $(Template_Include_Path) > gccTemplate/src/IPhreeqcMMS.d

# Compiles file src/IPhreeqcMMSLib.cpp for the Template configuration...
-include gccTemplate/src/IPhreeqcMMSLib.d
gccTemplate/src/IPhreeqcMMSLib.o: src/IPhreeqcMMSLib.cpp
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -c src/IPhreeqcMMSLib.cpp $(Template_Include_Path) -o gccTemplate/src/IPhreeqcMMSLib.o
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -MM src/IPhreeqcMMSLib.cpp $(Template_Include_Path) > gccTemplate/src/IPhreeqcMMSLib.d

# Compiles file src/stdcall.cpp for the Template configuration...
-include gccTemplate/src/stdcall.d
gccTemplate/src/stdcall.o: src/stdcall.cpp
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -c src/stdcall.cpp $(Template_Include_Path) -o gccTemplate/src/stdcall.o
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -MM src/stdcall.cpp $(Template_Include_Path) > gccTemplate/src/stdcall.d

# Builds the Template configuration...
.PHONY: Template
Template: create_folders x64/gccTemplate/src/cdecl.o x64/gccTemplate/src/fortran.o x64/gccTemplate/src/IPhreeqcMMS.o x64/gccTemplate/src/IPhreeqcMMSLib.o x64/gccTemplate/src/stdcall.o 
	g++ x64/gccTemplate/src/cdecl.o x64/gccTemplate/src/fortran.o x64/gccTemplate/src/IPhreeqcMMS.o x64/gccTemplate/src/IPhreeqcMMSLib.o x64/gccTemplate/src/stdcall.o  $(Template_Library_Path) $(Template_Libraries) -Wl,-rpath,./ -o ../x64/gccTemplate/IPhreeqcMMS.exe

# Compiles file src/cdecl.cpp for the Template configuration...
-include x64/gccTemplate/src/cdecl.d
x64/gccTemplate/src/cdecl.o: src/cdecl.cpp
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -c src/cdecl.cpp $(Template_Include_Path) -o x64/gccTemplate/src/cdecl.o
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -MM src/cdecl.cpp $(Template_Include_Path) > x64/gccTemplate/src/cdecl.d

# Compiles file src/fortran.cpp for the Template configuration...
-include x64/gccTemplate/src/fortran.d
x64/gccTemplate/src/fortran.o: src/fortran.cpp
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -c src/fortran.cpp $(Template_Include_Path) -o x64/gccTemplate/src/fortran.o
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -MM src/fortran.cpp $(Template_Include_Path) > x64/gccTemplate/src/fortran.d

# Compiles file src/IPhreeqcMMS.cpp for the Template configuration...
-include x64/gccTemplate/src/IPhreeqcMMS.d
x64/gccTemplate/src/IPhreeqcMMS.o: src/IPhreeqcMMS.cpp
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -c src/IPhreeqcMMS.cpp $(Template_Include_Path) -o x64/gccTemplate/src/IPhreeqcMMS.o
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -MM src/IPhreeqcMMS.cpp $(Template_Include_Path) > x64/gccTemplate/src/IPhreeqcMMS.d

# Compiles file src/IPhreeqcMMSLib.cpp for the Template configuration...
-include x64/gccTemplate/src/IPhreeqcMMSLib.d
x64/gccTemplate/src/IPhreeqcMMSLib.o: src/IPhreeqcMMSLib.cpp
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -c src/IPhreeqcMMSLib.cpp $(Template_Include_Path) -o x64/gccTemplate/src/IPhreeqcMMSLib.o
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -MM src/IPhreeqcMMSLib.cpp $(Template_Include_Path) > x64/gccTemplate/src/IPhreeqcMMSLib.d

# Compiles file src/stdcall.cpp for the Template configuration...
-include x64/gccTemplate/src/stdcall.d
x64/gccTemplate/src/stdcall.o: src/stdcall.cpp
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -c src/stdcall.cpp $(Template_Include_Path) -o x64/gccTemplate/src/stdcall.o
	$(CPP_COMPILER) $(Template_Preprocessor_Definitions) $(Template_Compiler_Flags) -MM src/stdcall.cpp $(Template_Include_Path) > x64/gccTemplate/src/stdcall.d

# Creates the intermediate and output folders for each configuration...
.PHONY: create_folders
create_folders:
	mkdir -p gccDebug/src
	mkdir -p ../gccDebug
	mkdir -p x64/gccDebug/src
	mkdir -p ../x64/gccDebug
	mkdir -p gccRelease/src
	mkdir -p ../gccRelease
	mkdir -p x64/gccRelease/src
	mkdir -p ../x64/gccRelease
	mkdir -p gccTemplate/src
	mkdir -p ../gccTemplate
	mkdir -p x64/gccTemplate/src
	mkdir -p ../x64/gccTemplate

# Cleans intermediate and output files (objects, libraries, executables)...
.PHONY: clean
clean:
	rm -f gccDebug/*.o
	rm -f gccDebug/*.d
	rm -f ../gccDebug/*.a
	rm -f ../gccDebug/*.so
	rm -f ../gccDebug/*.dll
	rm -f ../gccDebug/*.exe
	rm -f x64/gccDebug/*.o
	rm -f x64/gccDebug/*.d
	rm -f ../x64/gccDebug/*.a
	rm -f ../x64/gccDebug/*.so
	rm -f ../x64/gccDebug/*.dll
	rm -f ../x64/gccDebug/*.exe
	rm -f gccRelease/*.o
	rm -f gccRelease/*.d
	rm -f ../gccRelease/*.a
	rm -f ../gccRelease/*.so
	rm -f ../gccRelease/*.dll
	rm -f ../gccRelease/*.exe
	rm -f x64/gccRelease/*.o
	rm -f x64/gccRelease/*.d
	rm -f ../x64/gccRelease/*.a
	rm -f ../x64/gccRelease/*.so
	rm -f ../x64/gccRelease/*.dll
	rm -f ../x64/gccRelease/*.exe
	rm -f gccTemplate/*.o
	rm -f gccTemplate/*.d
	rm -f ../gccTemplate/*.a
	rm -f ../gccTemplate/*.so
	rm -f ../gccTemplate/*.dll
	rm -f ../gccTemplate/*.exe
	rm -f x64/gccTemplate/*.o
	rm -f x64/gccTemplate/*.d
	rm -f ../x64/gccTemplate/*.a
	rm -f ../x64/gccTemplate/*.so
	rm -f ../x64/gccTemplate/*.dll
	rm -f ../x64/gccTemplate/*.exe

