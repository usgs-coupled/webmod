# Compiler flags...
CPP_COMPILER = g++
C_COMPILER = gcc

# Include paths...
Debug_Include_Path=-I"src" -I"IPhreeqc/src" -I"IPhreeqc/src/phreeqcpp" 
Release_Include_Path=-I"src" -I"IPhreeqc/src" -I"IPhreeqc/src/phreeqcpp" 

# Library paths...
Debug_Library_Path=
Release_Library_Path=

# Additional libraries...
Debug_Libraries=
Release_Libraries=

# Preprocessor definitions...
Debug_Preprocessor_Definitions=-D _DEBUG -D GCC_BUILD -D _LIB -D SWIG_SHARED_OBJ -D _CRT_SECURE_NO_DEPRECATE 
Release_Preprocessor_Definitions=-D NDEBUG -D GCC_BUILD -D _LIB -D SWIG_SHARED_OBJ -D _CRT_SECURE_NO_DEPRECATE 

# Implictly linked object files...
Debug_Implicitly_Linked_Objects=
Release_Implicitly_Linked_Objects=

# Compiler flags...
Debug_Compiler_Flags=-O0 
Release_Compiler_Flags=-O2 

# Builds all configurations for this project...
.PHONY: build_all_configurations
build_all_configurations: Debug Release 

# Builds the Debug configuration...
.PHONY: Debug
Debug: create_folders gccDebug/src/cdecl.o gccDebug/src/fortran.o gccDebug/src/IPhreeqcMMS.o gccDebug/src/IPhreeqcMMSLib.o gccDebug/src/stdcall.o 
	ar rcs ../gccDebug/libIPhreeqcMMS.a gccDebug/src/cdecl.o gccDebug/src/fortran.o gccDebug/src/IPhreeqcMMS.o gccDebug/src/IPhreeqcMMSLib.o gccDebug/src/stdcall.o  $(Debug_Implicitly_Linked_Objects)


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


# Builds the Release configuration...
.PHONY: Release
Release: create_folders gccRelease/src/cdecl.o gccRelease/src/fortran.o gccRelease/src/IPhreeqcMMS.o gccRelease/src/IPhreeqcMMSLib.o gccRelease/src/stdcall.o 
	ar rcs ../gccRelease/libIPhreeqcMMS.a gccRelease/src/cdecl.o gccRelease/src/fortran.o gccRelease/src/IPhreeqcMMS.o gccRelease/src/IPhreeqcMMSLib.o gccRelease/src/stdcall.o  $(Release_Implicitly_Linked_Objects)

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


# Creates the intermediate and output folders for each configuration...
.PHONY: create_folders
create_folders:
	mkdir -p gccDebug/src
	mkdir -p ../gccDebug
	mkdir -p gccRelease/src
	mkdir -p ../gccRelease

# Cleans intermediate and output files (objects, libraries, executables)...
.PHONY: clean
clean:
	rm -f gccDebug/*.o
	rm -f gccDebug/*.d
	rm -f ../gccDebug/*.a
	rm -f ../gccDebug/*.so
	rm -f ../gccDebug/*.dll
	rm -f ../gccDebug/*.exe
	rm -f gccRelease/*.o
	rm -f gccRelease/*.d
	rm -f ../gccRelease/*.a
	rm -f ../gccRelease/*.so
	rm -f ../gccRelease/*.dll
	rm -f ../gccRelease/*.exe

