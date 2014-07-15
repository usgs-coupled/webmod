# Compiler flags...
CPP_COMPILER = g++
C_COMPILER = gcc
F_COMPILER = /opt/intel/bin/ifort

# Include paths...
Debug_Include_Path=-I../IPhreeqcMMS/include -I../IPhreeqcMMS/IPhreeqc/src
Release_Include_Path=-I../IPhreeqcMMS/include -I../IPhreeqcMMS/IPhreeqc/src

# Library paths...
Debug_Library_Path=
Release_Library_Path=

# Additional libraries...
Debug_Libraries=
Release_Libraries=

# Preprocessor definitions...
Debug_Preprocessor_Definitions=-D GCC_BUILDXX -D _DEBUG -D _LIB -D _CRT_SECURE_NO_DEPRECATE -D _CRT_NONSTDC_NO_DEPRECATE 
Release_Preprocessor_Definitions=-D GCC_BUILDXX -D NDEBUG -D _LIB -D _CRT_SECURE_NO_DEPRECATE -D _CRT_NONSTDC_NO_DEPRECATE 

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
Debug: 

gccDebug/src/basin_topg.o: src/basin_topg.f src/fmodules.inc                             #       MODULE WEBMOD_BASIN
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/basin_topg.f \
		$(Debug_Include_Path) -o gccDebug/src/basin_topg.o

gccDebug/src/call_modules_webmod.o: src/call_modules_webmod.f src/fmodules.inc           #       MODULE WEB_MODULE
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/call_modules_webmod.f \
		$(Debug_Include_Path) -o gccDebug/src/call_modules_webmod.o

gccDebug/src/ccsolrad_web.o: src/ccsolrad_web.f src/fmodules.inc                         #       MODULE WEBMOD_CSOL
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/ccsolrad_web.f \
		$(Debug_Include_Path) -o gccDebug/src/ccsolrad_web.o

gccDebug/src/elapsed_time_function.o: src/elapsed_time_function.f90
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/elapsed_time_function.f90 \
		$(Debug_Include_Path) -o gccDebug/src/elapsed_time_function.o

gccDebug/src/intcp_prms.o: src/intcp_prms.f src/fmodules.inc                             #       MODULE WEBMOD_INTCP
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/intcp_prms.f \
		$(Debug_Include_Path) -o gccDebug/src/intcp_prms.o

gccDebug/src/io.o: src/io.f src/fmodules.inc                                             #       MODULE WEBMOD_IO
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/io.f \
		$(Debug_Include_Path) -o gccDebug/src/io.o

gccDebug/src/irrig_web.o: src/irrig_web.f src/fmodules.inc                               #       MODULE WEBMOD_IRRIG
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/irrig_web.f \
		$(Debug_Include_Path) -o gccDebug/src/irrig_web.o

gccDebug/src/nwsmelt_topg.o: src/nwsmelt_topg.f src/fmodules.inc \
	gccDebug/src/obs_webmod.o                                                #       MODULE WEBMOD_SNOW
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/nwsmelt_topg.f \
		$(Debug_Include_Path) -o gccDebug/src/nwsmelt_topg.o

gccDebug/src/obs_chem.o: src/obs_chem.f  src/fmodules.inc                                #       MODULE WEBMOD_OBSCHEM
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/obs_chem.f \
		$(Debug_Include_Path) -o gccDebug/src/obs_chem.o

gccDebug/src/obs_webmod.o: src/obs_webmod.f src/fmodules.inc                             #       MODULE WEBMOD_OBSHYD
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/obs_webmod.f \
		$(Debug_Include_Path) -o gccDebug/src/obs_webmod.o

gccDebug/src/obs_webx.o: src/obs_webx.f src/fmodules.inc                             
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/obs_webx.f \
		$(Debug_Include_Path) -o gccDebug/src/obs_webx.o

gccDebug/src/phreeq_mms.o: src/phreeq_mms.f90 src/fmodules.inc ../IPhreeqcMMS/IPhreeqc/src/IPhreeqc.f90.inc \
	../IPhreeqcMMS/include/mms_phreeqc.f90.inc \
	gccDebug/src/io.o gccDebug/src/obs_chem.o \
	gccDebug/src/obs_webmod.o \
	gccDebug/src/irrig_web.o gccDebug/src/potet_hamon_prms.o \
	gccDebug/src/temp_1sta_prms.o gccDebug/src/webmod_res.o \
	gccDebug/src/topmod_chem.o                                                            #      MODULE WEBMOD_PHREEQ_MMS
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/phreeq_mms.f90 \
		$(Debug_Include_Path) -o gccDebug/src/phreeq_mms.o

gccDebug/src/potet_hamon_prms.o: src/potet_hamon_prms.f src/fmodules.inc                              #       MODULE WEBMOD_POTET
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/potet_hamon_prms.f \
		$(Debug_Include_Path) -o gccDebug/src/potet_hamon_prms.o

gccDebug/src/precip_web.o: src/precip_web.f src/fmodules.inc                                          #       MODULE WEBMOD_PRECIP
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/precip_web.f \
		$(Debug_Include_Path) -o gccDebug/src/precip_web.o

gccDebug/src/route_clark.o: src/route_clark.f src/fmodules.inc gccDebug/src/io.o                      #       MODULE WEBMOD_ROUTE
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/route_clark.f \
		$(Debug_Include_Path) -o gccDebug/src/route_clark.o

gccDebug/src/soltab_prms.o: src/soltab_prms.f src/fmodules.inc                                        #       MODULE WEBMOD_SOLTAB
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/soltab_prms.f \
		$(Debug_Include_Path) -o gccDebug/src/soltab_prms.o

gccDebug/src/temp_1sta_prms.o: src/temp_1sta_prms.f src/fmodules.inc                                  #       MODULE WEBMOD_TEMP1STA
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/temp_1sta_prms.f \
		$(Debug_Include_Path) -o gccDebug/src/temp_1sta_prms.o

gccDebug/src/top2clark.o: src/top2clark.f src/fmodules.inc gccDebug/src/io.o                          #       MODULE WEBMOD_CLARK
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/top2clark.f \
		$(Debug_Include_Path) -o gccDebug/src/top2clark.o

gccDebug/src/topmod_chem.o: src/topmod_chem.f src/fmodules.inc gccDebug/src/io.o                      #       MODULE WEBMOD_TOPMOD
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/topmod_chem.f \
		$(Debug_Include_Path) -o gccDebug/src/topmod_chem.o

gccDebug/src/webmod_res.o: src/webmod_res.f src/fmodules.inc gccDebug/src/topmod_chem.o               #       MODULE WEBMOD_RESMOD
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/webmod_res.f \
		$(Debug_Include_Path) -o gccDebug/src/webmod_res.o

gccDebug/src/web_sum.o: src/web_sum.F src/fmodules.inc gccDebug/src/webmod_res.o                      #       MODULE WEBMOD_SUM
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/web_sum.F \
		$(Debug_Include_Path) -o gccDebug/src/web_sum.o

gccDebug/src/xyz_dist_web.o: src/xyz_dist_web.f src/fmodules.inc 
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c src/xyz_dist_web.f \
		$(Debug_Include_Path) -o gccDebug/src/xyz_dist_web.o

gccDebug/src/phr_mix.o: ../IPhreeqcMMS/src/phr_mix.f90 gccDebug/src/io.o gccDebug/src/phreeq_mms.o
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c ../IPhreeqcMMS/src/phr_mix.f90 \
		$(Debug_Include_Path) -o gccDebug/src/phr_mix.o

gccDebug/src/phr_multicopy.o: ../IPhreeqcMMS/src/phr_multicopy.f90
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c ../IPhreeqcMMS/src/phr_multicopy.f90 \
		$(Debug_Include_Path) -o gccDebug/src/phr_multicopy.o

gccDebug/src/phr_precip.o: ../IPhreeqcMMS/src/phr_precip.f90
	$(F_COMPILER) $(Debug_Preprocessor_Definitions) $(Debug_Compiler_Flags) -c ../IPhreeqcMMS/src/phr_precip.f90 \
		$(Debug_Include_Path) -o gccDebug/src/phr_precip.o

# Builds the Debug configuration...
.PHONY: Debug
Debug: create_folders  \
gccDebug/src/basin_topg.o \
gccDebug/src/call_modules_webmod.o \
gccDebug/src/ccsolrad_web.o \
gccDebug/src/elapsed_time_function.o \
gccDebug/src/intcp_prms.o \
gccDebug/src/io.o \
gccDebug/src/irrig_web.o \
gccDebug/src/nwsmelt_topg.o \
gccDebug/src/obs_chem.o \
gccDebug/src/obs_webmod.o \
gccDebug/src/obs_webx.o \
gccDebug/src/phreeq_mms.o \
gccDebug/src/potet_hamon_prms.o \
gccDebug/src/precip_web.o \
gccDebug/src/route_clark.o \
gccDebug/src/soltab_prms.o \
gccDebug/src/temp_1sta_prms.o \
gccDebug/src/top2clark.o \
gccDebug/src/topmod_chem.o \
gccDebug/src/webmod_res.o \
gccDebug/src/web_sum.o \
gccDebug/src/xyz_dist_web.o \
gccDebug/src/phr_mix.o \
gccDebug/src/phr_multicopy.o \
gccDebug/src/phr_precip.o


# Builds the Release configuration...
.PHONY: Release
Release: 

gccRelease/src/basin_topg.o: src/basin_topg.f src/fmodules.inc                             #       MODULE WEBMOD_BASIN
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/basin_topg.f \
		$(Release_Include_Path) -o gccRelease/src/basin_topg.o

gccRelease/src/call_modules_webmod.o: src/call_modules_webmod.f src/fmodules.inc           #       MODULE WEB_MODULE
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/call_modules_webmod.f \
		$(Release_Include_Path) -o gccRelease/src/call_modules_webmod.o

gccRelease/src/ccsolrad_web.o: src/ccsolrad_web.f src/fmodules.inc                         #       MODULE WEBMOD_CSOL
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/ccsolrad_web.f \
		$(Release_Include_Path) -o gccRelease/src/ccsolrad_web.o

gccRelease/src/elapsed_time_function.o: src/elapsed_time_function.f90
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/elapsed_time_function.f90 \
		$(Release_Include_Path) -o gccRelease/src/elapsed_time_function.o

gccRelease/src/intcp_prms.o: src/intcp_prms.f src/fmodules.inc                             #       MODULE WEBMOD_INTCP
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/intcp_prms.f \
		$(Release_Include_Path) -o gccRelease/src/intcp_prms.o

gccRelease/src/io.o: src/io.f src/fmodules.inc                                             #       MODULE WEBMOD_IO
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/io.f \
		$(Release_Include_Path) -o gccRelease/src/io.o

gccRelease/src/irrig_web.o: src/irrig_web.f src/fmodules.inc                               #       MODULE WEBMOD_IRRIG
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/irrig_web.f \
		$(Release_Include_Path) -o gccRelease/src/irrig_web.o

gccRelease/src/nwsmelt_topg.o: src/nwsmelt_topg.f src/fmodules.inc \
	gccRelease/src/obs_webmod.o                                                #       MODULE WEBMOD_SNOW
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/nwsmelt_topg.f \
		$(Release_Include_Path) -o gccRelease/src/nwsmelt_topg.o

gccRelease/src/obs_chem.o: src/obs_chem.f  src/fmodules.inc                                #       MODULE WEBMOD_OBSCHEM
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/obs_chem.f \
		$(Release_Include_Path) -o gccRelease/src/obs_chem.o

gccRelease/src/obs_webmod.o: src/obs_webmod.f src/fmodules.inc                             #       MODULE WEBMOD_OBSHYD
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/obs_webmod.f \
		$(Release_Include_Path) -o gccRelease/src/obs_webmod.o

gccRelease/src/obs_webx.o: src/obs_webx.f src/fmodules.inc                             
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/obs_webx.f \
		$(Release_Include_Path) -o gccRelease/src/obs_webx.o

gccRelease/src/phreeq_mms.o: src/phreeq_mms.f90 src/fmodules.inc ../IPhreeqcMMS/IPhreeqc/src/IPhreeqc.f90.inc \
	../IPhreeqcMMS/include/mms_phreeqc.f90.inc \
	gccRelease/src/io.o gccRelease/src/obs_chem.o \
	gccRelease/src/obs_webmod.o \
	gccRelease/src/irrig_web.o gccRelease/src/potet_hamon_prms.o \
	gccRelease/src/temp_1sta_prms.o gccRelease/src/webmod_res.o \
	gccRelease/src/topmod_chem.o                                                            #      MODULE WEBMOD_PHREEQ_MMS
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/phreeq_mms.f90 \
		$(Release_Include_Path) -o gccRelease/src/phreeq_mms.o

gccRelease/src/potet_hamon_prms.o: src/potet_hamon_prms.f src/fmodules.inc                              #       MODULE WEBMOD_POTET
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/potet_hamon_prms.f \
		$(Release_Include_Path) -o gccRelease/src/potet_hamon_prms.o

gccRelease/src/precip_web.o: src/precip_web.f src/fmodules.inc                                          #       MODULE WEBMOD_PRECIP
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/precip_web.f \
		$(Release_Include_Path) -o gccRelease/src/precip_web.o

gccRelease/src/route_clark.o: src/route_clark.f src/fmodules.inc gccRelease/src/io.o                      #       MODULE WEBMOD_ROUTE
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/route_clark.f \
		$(Release_Include_Path) -o gccRelease/src/route_clark.o

gccRelease/src/soltab_prms.o: src/soltab_prms.f src/fmodules.inc                                        #       MODULE WEBMOD_SOLTAB
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/soltab_prms.f \
		$(Release_Include_Path) -o gccRelease/src/soltab_prms.o

gccRelease/src/temp_1sta_prms.o: src/temp_1sta_prms.f src/fmodules.inc                                  #       MODULE WEBMOD_TEMP1STA
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/temp_1sta_prms.f \
		$(Release_Include_Path) -o gccRelease/src/temp_1sta_prms.o

gccRelease/src/top2clark.o: src/top2clark.f src/fmodules.inc gccRelease/src/io.o                          #       MODULE WEBMOD_CLARK
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/top2clark.f \
		$(Release_Include_Path) -o gccRelease/src/top2clark.o

gccRelease/src/topmod_chem.o: src/topmod_chem.f src/fmodules.inc gccRelease/src/io.o                      #       MODULE WEBMOD_TOPMOD
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/topmod_chem.f \
		$(Release_Include_Path) -o gccRelease/src/topmod_chem.o

gccRelease/src/webmod_res.o: src/webmod_res.f src/fmodules.inc gccRelease/src/topmod_chem.o               #       MODULE WEBMOD_RESMOD
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/webmod_res.f \
		$(Release_Include_Path) -o gccRelease/src/webmod_res.o

gccRelease/src/web_sum.o: src/web_sum.F src/fmodules.inc gccRelease/src/webmod_res.o                      #       MODULE WEBMOD_SUM
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/web_sum.F \
		$(Release_Include_Path) -o gccRelease/src/web_sum.o

gccRelease/src/xyz_dist_web.o: src/xyz_dist_web.f src/fmodules.inc 
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c src/xyz_dist_web.f \
		$(Release_Include_Path) -o gccRelease/src/xyz_dist_web.o

gccRelease/src/phr_mix.o: ../IPhreeqcMMS/src/phr_mix.f90 gccRelease/src/io.o gccRelease/src/phreeq_mms.o
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c ../IPhreeqcMMS/src/phr_mix.f90 \
		$(Release_Include_Path) -o gccRelease/src/phr_mix.o

gccRelease/src/phr_multicopy.o: ../IPhreeqcMMS/src/phr_multicopy.f90
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c ../IPhreeqcMMS/src/phr_multicopy.f90 \
		$(Release_Include_Path) -o gccRelease/src/phr_multicopy.o

gccRelease/src/phr_precip.o: ../IPhreeqcMMS/src/phr_precip.f90
	$(F_COMPILER) $(Release_Preprocessor_Definitions) $(Release_Compiler_Flags) -c ../IPhreeqcMMS/src/phr_precip.f90 \
		$(Release_Include_Path) -o gccRelease/src/phr_precip.o

# Builds the Release configuration...
.PHONY: Release
Release: create_folders  \
gccRelease/src/basin_topg.o \
gccRelease/src/call_modules_webmod.o \
gccRelease/src/ccsolrad_web.o \
gccRelease/src/elapsed_time_function.o \
gccRelease/src/intcp_prms.o \
gccRelease/src/io.o \
gccRelease/src/irrig_web.o \
gccRelease/src/nwsmelt_topg.o \
gccRelease/src/obs_chem.o \
gccRelease/src/obs_webmod.o \
gccRelease/src/obs_webx.o \
gccRelease/src/phreeq_mms.o \
gccRelease/src/potet_hamon_prms.o \
gccRelease/src/precip_web.o \
gccRelease/src/route_clark.o \
gccRelease/src/soltab_prms.o \
gccRelease/src/temp_1sta_prms.o \
gccRelease/src/top2clark.o \
gccRelease/src/topmod_chem.o \
gccRelease/src/webmod_res.o \
gccRelease/src/web_sum.o \
gccRelease/src/xyz_dist_web.o \
gccRelease/src/phr_mix.o \
gccRelease/src/phr_multicopy.o \
gccRelease/src/phr_precip.o

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
	rm -f gccDebug/src/*.o
	rm -f gccDebug/src/*.d
	rm -f ../gccDebug/*.a
	rm -f ../gccDebug/*.so
	rm -f ../gccDebug/*.dll
	rm -f ../gccDebug/*.exe
	rm -f gccRelease/src/*.o
	rm -f gccRelease/src/*.d
	rm -f ../gccRelease/*.a
	rm -f ../gccRelease/*.so
	rm -f ../gccRelease/*.dll
	rm -f ../gccRelease/*.exe
