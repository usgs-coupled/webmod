EXTRA_DIST=\
	CMakeLists.txt\
	README_WEBMOD.TXT\
	doc/CMakeLists.txt\
	mmf_c/CMakeLists.txt\
	IPhreeqcMMS/CMakeLists.txt\
	IPhreeqcMMS/IPhreeqc/CMakeLists.txt\
	lib/CMakeLists.txt\
	Tests/CMakeLists.txt\
	Tests/Andrews/CMakeLists.txt\
	Tests/Andrews/control/andcrk.control.in\
	Tests/Andrews/input/andcrk.params\
	Tests/Andrews/input/andcrk.dat\
	Tests/Andrews/input/phreeq_lut\
	Tests/Andrews/input/andcrk.dat.chemdat\
	Tests/Andrews/input/webmod.pqi\
	Tests/Andrews/input/webmod.chem.dat\
	Tests/Andrews/input/phreeqc_web_lite.dat\
	webmod.d/CMakeLists.txt


# $(ANDREWS)\
# $(ANDREWS_CONTROL)\
# $(ANDREWS_INPUT)\
# $(ANDREWS_OUTPUT)\
# $(ANDREWS_TUTORIAL)\
# $(ANDREWS_TUTORIAL_CONTROL)\
# $(ANDREWS_TUTORIAL_INPUT)\
# $(ANDREWS_TUTORIAL_OUTPUT)\
# $(DR2)\
# $(DR2_CONTROL)\
# $(DR2_INPUT)\
# $(DR2_OUTPUT)


SUBDIRS = lib

AM_CPPFLAGS=-I$(srcdir)/IPhreeqcMMS/src -I$(srcdir)/IPhreeqcMMS/IPhreeqc/src -I$(srcdir)/IPhreeqcMMS/IPhreeqc/src/phreeqcpp -I$(srcdir)/IPhreeqcMMS/IPhreeqc/src/phreeqcpp/common -I$(srcdir)/IPhreeqcMMS/IPhreeqc/src/phreeqcpp/PhreeqcKeywords

AM_FCFLAGS=-I$(srcdir)/IPhreeqcMMS/include

FCDEFS=

# doc
dist_doc_DATA = $(DOCS)

DOCS=\
	Distribution/doc/Webb_WEBMOD_Body_2016Aug19.docx\
	Distribution/doc/Webb_WEBMOD_Tables_2016August19.docm\
	Distribution/doc/WEBMOD_figs_1_12_2016August19.pdf\
	Distribution/doc/WEBMOD_figs_13_52_2016August19.pdf

# Andrews
doc_projects_Andrewsdir = $(docdir)/projects/Andrews

dist_doc_projects_Andrews_DATA = $(ANDREWS)

ANDREWS=\
	Distribution/projects/Andrews/Andrews.xlsm\
	Distribution/projects/Andrews/webmod.bat\
	Distribution/projects/Andrews/webmod_gui.bat\
	Distribution/projects/Andrews/webmod_paramtool.bat\
	Distribution/projects/Andrews/webmod_print.bat

# Andrews/control
doc_projects_Andrews_controldir = $(docdir)/projects/Andrews/control

dist_doc_projects_Andrews_control_DATA = $(ANDREWS_CONTROL)

ANDREWS_CONTROL=\
	Distribution/projects/Andrews/control/webmod.control

# Andrews/input
doc_projects_Andrews_inputdir = $(docdir)/projects/Andrews/input

dist_doc_projects_Andrews_input_DATA = $(ANDREWS_INPUT)

ANDREWS_INPUT=\
	Distribution/projects/Andrews/input/phreeqc_web_lite.dat\
	Distribution/projects/Andrews/input/phreeq_lut\
	Distribution/projects/Andrews/input/webmod.chem.dat\
	Distribution/projects/Andrews/input/webmod.hydro.dat\
	Distribution/projects/Andrews/input/webmod.params\
	Distribution/projects/Andrews/input/webmod.pqi

# Andrews/output
doc_projects_Andrews_outputdir = $(docdir)/projects/Andrews/output

dist_doc_projects_Andrews_output_DATA = $(ANDREWS_OUTPUT)

ANDREWS_OUTPUT=\
	Distribution/projects/Andrews/output/webmod.chem.out\
	Distribution/projects/Andrews/output/webmod.hydro.out\
	Distribution/projects/Andrews/output/webmod.statvar\
	Distribution/projects/Andrews/output/webmod.topout


# Andrews_tutorial
doc_projects_Andrews_tutorialdir = $(docdir)/projects/Andrews_tutorial

dist_doc_projects_Andrews_tutorial_DATA = $(ANDREWS_TUTORIAL)

ANDREWS_TUTORIAL=\
	Distribution/projects/Andrews_tutorial/webmod.bat\
	Distribution/projects/Andrews_tutorial/webmod_gui.bat\
	Distribution/projects/Andrews_tutorial/webmod_paramtool.bat\
	Distribution/projects/Andrews_tutorial/webmod_print.bat

# Andrews_tutorial/control
doc_projects_Andrews_tutorial_controldir = $(docdir)/projects/Andrews_tutorial/control

dist_doc_projects_Andrews_tutorial_control_DATA = $(ANDREWS_TUTORIAL_CONTROL)

ANDREWS_TUTORIAL_CONTROL=\
	Distribution/projects/Andrews_tutorial/control/webmod.control

# Andrews_tutorial/input
doc_projects_Andrews_tutorial_inputdir = $(docdir)/projects/Andrews_tutorial/input

dist_doc_projects_Andrews_tutorial_input_DATA = $(ANDREWS_TUTORIAL_INPUT)

ANDREWS_TUTORIAL_INPUT=\
	Distribution/projects/Andrews_tutorial/input/phreeqc_web_lite.dat\
	Distribution/projects/Andrews_tutorial/input/phreeq_lut\
	Distribution/projects/Andrews_tutorial/input/webmod.chem.dat\
	Distribution/projects/Andrews_tutorial/input/webmod.hydro.dat\
	Distribution/projects/Andrews_tutorial/input/webmod.params\
	Distribution/projects/Andrews_tutorial/input/webmod.pqi

# Andrews_tutorial/output
doc_projects_Andrews_tutorial_outputdir = $(docdir)/projects/Andrews_tutorial/output

dist_doc_projects_Andrews_tutorial_output_DATA = $(ANDREWS_TUTORIAL_OUTPUT)

ANDREWS_TUTORIAL_OUTPUT=\
	Distribution/projects/Andrews_tutorial/output/webmod.chem.out\
	Distribution/projects/Andrews_tutorial/output/webmod.hydro.out\
	Distribution/projects/Andrews_tutorial/output/webmod.statvar\
	Distribution/projects/Andrews_tutorial/output/webmod.topout


# dr2
doc_projects_dr2dir = $(docdir)/projects/dr2

dist_doc_projects_dr2_DATA = $(DR2)

DR2=\
	Distribution/projects/dr2/dr2.xlsm\
	Distribution/projects/dr2/webmod.bat\
	Distribution/projects/dr2/webmod_gui.bat\
	Distribution/projects/dr2/webmod_paramtool.bat\
	Distribution/projects/dr2/webmod_print.bat


# dr2/control
doc_projects_dr2_controldir = $(docdir)/projects/dr2/control

dist_doc_projects_dr2_control_DATA = $(DR2_CONTROL)

DR2_CONTROL=\
	Distribution/projects/dr2/control/webmod.control

# dr2/input
doc_projects_dr2_inputdir = $(docdir)/projects/dr2/input

dist_doc_projects_dr2_input_DATA = $(DR2_INPUT)

DR2_INPUT=\
	Distribution/projects/dr2/input/phreeqc_web_lite.dat\
	Distribution/projects/dr2/input/phreeq_lut\
	Distribution/projects/dr2/input/webmod.hydro.dat\
	Distribution/projects/dr2/input/webmod.params\
	Distribution/projects/dr2/input/webmod.pqi

# dr2/output
doc_projects_dr2_outputdir = $(docdir)/projects/dr2/output

dist_doc_projects_dr2_output_DATA = $(DR2_OUTPUT)

DR2_OUTPUT=\
	Distribution/projects/dr2/output/webmod.chem.out\
	Distribution/projects/dr2/output/webmod.hydro.out\
	Distribution/projects/dr2/output/webmod.statvar\
	Distribution/projects/dr2/output/webmod.topout

# bin
bin_PROGRAMS = webmod

webmod_SOURCES=\
	.depend\
	IPhreeqcMMS/src/cdecl.cpp\
	IPhreeqcMMS/src/fortran.cpp\
	IPhreeqcMMS/src/fortran.h\
	IPhreeqcMMS/src/IPhreeqcMMS.cpp\
	IPhreeqcMMS/src/IPhreeqcMMS.h\
	IPhreeqcMMS/src/IPhreeqcMMS.hpp\
	IPhreeqcMMS/src/IPhreeqcMMSLib.cpp\
	IPhreeqcMMS/src/IPhreeqcMMSLib.h\
	IPhreeqcMMS/include/mms_phreeqc.f90.inc\
	IPhreeqcMMS/src/phr_mix.F90\
	IPhreeqcMMS/src/phr_multicopy.f90\
	IPhreeqcMMS/src/phr_precip.f90\
	IPhreeqcMMS/src/stdcall.cpp\
	IPhreeqcMMS/IPhreeqc/src/CSelectedOutput.cpp\
	IPhreeqcMMS/IPhreeqc/src/CSelectedOutput.hxx\
	IPhreeqcMMS/IPhreeqc/src/CVar.hxx\
	IPhreeqcMMS/IPhreeqc/src/Debug.h\
	IPhreeqcMMS/IPhreeqc/src/ErrorReporter.hxx\
	IPhreeqcMMS/IPhreeqc/src/fwrap.cpp\
	IPhreeqcMMS/IPhreeqc/src/fwrap.h\
	IPhreeqcMMS/IPhreeqc/src/IPhreeqc.cpp\
	IPhreeqcMMS/IPhreeqc/src/IPhreeqc.f.inc\
	IPhreeqcMMS/IPhreeqc/src/IPhreeqc.f90.inc\
	IPhreeqcMMS/IPhreeqc/src/IPhreeqc.h\
	IPhreeqcMMS/IPhreeqc/src/IPhreeqc.hpp\
	IPhreeqcMMS/IPhreeqc/src/IPhreeqc_interface.F90\
	IPhreeqcMMS/IPhreeqc/src/IPhreeqc_interface_F.cpp\
	IPhreeqcMMS/IPhreeqc/src/IPhreeqc_interface_F.h\
	IPhreeqcMMS/IPhreeqc/src/IPhreeqcCallbacks.h\
	IPhreeqcMMS/IPhreeqc/src/IPhreeqcLib.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/advection.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/basicsubs.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/ChartHandler.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/ChartHandler.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/ChartObject.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/ChartObject.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/cl1.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/common/Parser.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/common/Parser.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/common/PHRQ_base.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/common/PHRQ_base.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/common/PHRQ_io.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/common/PHRQ_io.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/common/phrqtype.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/common/Utils.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/common/Utils.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/CurveObject.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/CurveObject.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/cvdense.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/cvdense.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/cvode.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/cvode.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/cxxKinetics.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/cxxKinetics.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/cxxMix.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/cxxMix.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/dense.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/dense.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Dictionary.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Dictionary.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/dumper.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/dumper.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Exchange.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Exchange.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/ExchComp.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/ExchComp.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/GasComp.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/GasComp.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/gases.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/GasPhase.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/GasPhase.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/global_structures.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/input.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/integrate.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/inverse.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/ISolution.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/ISolution.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/ISolutionComp.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/ISolutionComp.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/isotopes.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/kinetics.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/KineticsComp.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/KineticsComp.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/mainsubs.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/model.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/NA.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/NameDouble.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/NameDouble.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/NumKeyword.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/NumKeyword.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/nvector.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/nvector.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/nvector_serial.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/nvector_serial.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/parse.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/PBasic.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/PBasic.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/phqalloc.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/phqalloc.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Phreeqc.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Phreeqc.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/PhreeqcKeywords/Keywords.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/PhreeqcKeywords/Keywords.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/PHRQ_io_output.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/pitzer.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/pitzer_structures.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/PPassemblage.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/PPassemblage.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/PPassemblageComp.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/PPassemblageComp.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/prep.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Pressure.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Pressure.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/print.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Reaction.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Reaction.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/read.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/ReadClass.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/readtr.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/runner.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/runner.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SelectedOutput.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SelectedOutput.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Serializer.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Serializer.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/sit.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/smalldense.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/smalldense.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Solution.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Solution.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SolutionIsotope.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SolutionIsotope.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/spread.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SS.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SS.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SSassemblage.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SSassemblage.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SScomp.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SScomp.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/step.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/StorageBin.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/StorageBin.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/StorageBinList.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/StorageBinList.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/structures.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/sundialsmath.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/sundialsmath.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/sundialstypes.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Surface.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Surface.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SurfaceCharge.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SurfaceCharge.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SurfaceComp.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/SurfaceComp.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/System.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/System.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/tally.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Temperature.cxx\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Temperature.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/tidy.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/transport.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Use.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/Use.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/UserPunch.cpp\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/UserPunch.h\
	IPhreeqcMMS/IPhreeqc/src/phreeqcpp/utilities.cpp\
	IPhreeqcMMS/IPhreeqc/src/thread.h\
	IPhreeqcMMS/IPhreeqc/src/Var.c\
	IPhreeqcMMS/IPhreeqc/src/Var.h\
	IPhreeqcMMS/IPhreeqc/src/Version.h\
	webmod.d/src/basin_topg.F\
	webmod.d/src/call_modules_webmod.F\
	webmod.d/src/ccsolrad_web.F\
	webmod.d/src/defines.h\
	webmod.d/src/elapsed_time_function.F90\
	webmod.d/src/fmodules.inc\
	webmod.d/src/fmodules.f90.inc\
	webmod.d/src/intcp_prms.F\
	webmod.d/src/io.F\
	webmod.d/src/irrig_web.F\
	webmod.d/src/my_newunit.F\
	webmod.d/src/nwsmelt_topg.F\
	webmod.d/src/obs_chem.F\
	webmod.d/src/obs_webmod.F\
	webmod.d/src/phreeq_mms.F90\
	webmod.d/src/potet_hamon_prms.F\
	webmod.d/src/precip_web.F\
	webmod.d/src/route_clark.F\
	webmod.d/src/soltab_prms.F\
	webmod.d/src/temp_1sta_prms.F\
	webmod.d/src/top2clark.F\
	webmod.d/src/topmod_chem.F\
	webmod.d/src/web_sum.F\
	webmod.d/src/webmod_res.F\
	mmf_c/src/alloc_space.c\
	mmf_c/src/batch_run.c\
	mmf_c/src/batch_run_functions.c\
	mmf_c/src/build_lists.c\
	mmf_c/src/call_modules.c\
	mmf_c/src/call_setdims.c\
	mmf_c/src/check_vars.c\
	mmf_c/src/control_addr.c\
	mmf_c/src/control_array.c\
	mmf_c/src/control_var.c\
	mmf_c/src/create_vstats.c\
	mmf_c/src/decl_control.c\
	mmf_c/src/decldim.c\
	mmf_c/src/declparam.c\
	mmf_c/src/declvar.c\
	mmf_c/src/defs.h\
	mmf_c/src/dim_addr.c\
	mmf_c/src/dprint.c\
	mmf_c/src/free_vstats.c\
	mmf_c/src/get_elem_add.c\
	mmf_c/src/get_times.c\
	mmf_c/src/getdim.c\
	mmf_c/src/getdimname.c\
	mmf_c/src/getparam.c\
	mmf_c/src/getvar.c\
	mmf_c/src/globals.h\
	mmf_c/src/graph_single_run.c\
	mmf_c/src/julconvert.c\
	mmf_c/src/julday.c\
	mmf_c/src/load_param.c\
	mmf_c/src/mmf.c\
	mmf_c/src/mms.h\
	mmf_c/src/nodes.h\
	mmf_c/src/oprint.c\
	mmf_c/src/param_addr.c\
	mmf_c/src/parse_args.c\
	mmf_c/src/print_model_info.c\
	mmf_c/src/print_params.c\
	mmf_c/src/print_vars.c\
	mmf_c/src/protos.h\
	mmf_c/src/putvar.c\
	mmf_c/src/read_control.c\
	mmf_c/src/read_datainfo.c\
	mmf_c/src/read_line.c\
	mmf_c/src/read_params.c\
	mmf_c/src/read_vars.c\
	mmf_c/src/readvar.c\
	mmf_c/src/reset_dim.c\
	mmf_c/src/save_params.c\
	mmf_c/src/save_vars.c\
	mmf_c/src/setup_cont.c\
	mmf_c/src/sort_dims.c\
	mmf_c/src/sort_params.c\
	mmf_c/src/sort_vars.c\
	mmf_c/src/stats.c\
	mmf_c/src/str_to_vals.c\
	mmf_c/src/structs.h\
	mmf_c/src/timing.c\
	mmf_c/src/umalloc_etc.c\
	mmf_c/src/uprint.c\
	mmf_c/src/var_addr.c\
	mmf_c/src/write_vstats.c


webmod_LDADD = $(FCLIBS)

# We must override the default PPFCCOMPILE since ifort can't
# handle the FC_FUNC define correctly so $(DEFS) is replaced
# by $(FCDEFS)
#
PPFCCOMPILE = $(FC) $(FCDEFS) $(DEFAULT_INCLUDES) $(INCLUDES) \
	$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_FCFLAGS) $(FCFLAGS)

CLEANFILES=\
	iphreeqc.$(F90_MODEXT)\
	web_module.$(F90_MODEXT)\
	webmod_basin.$(F90_MODEXT)\
	webmod_clark.$(F90_MODEXT)\
	webmod_csol.$(F90_MODEXT)\
	webmod_intcp.$(F90_MODEXT)\
	webmod_io.$(F90_MODEXT)\
	webmod_irrig.$(F90_MODEXT)\
	webmod_obschem.$(F90_MODEXT)\
	webmod_obshyd.$(F90_MODEXT)\
	webmod_phreeq_mms.$(F90_MODEXT)\
	webmod_potet.$(F90_MODEXT)\
	webmod_precip.$(F90_MODEXT)\
	webmod_resmod.$(F90_MODEXT)\
	webmod_route.$(F90_MODEXT)\
	webmod_snow.$(F90_MODEXT)\
	webmod_soltab.$(F90_MODEXT)\
	webmod_sum.$(F90_MODEXT)\
	webmod_temp1sta.$(F90_MODEXT)\
	webmod_topmod.$(F90_MODEXT)

include .depend
