TESTS = Andrews

EXTRA_DIST = \
	$(TESTS) \
	CMakeLists.txt \
	input/andcrk.dat \
	input/andcrk.dat.chemdat \
	input/andcrk.params \
	input/phreeqc_web_lite.dat \
	input/phreeq_lut \
	input/webmod.chem.dat \
	input/webmod.pqi

TESTS_ENVIRONMENT = top_srcdir=$(top_srcdir) \
                    PATH=.$(PATH_SEPARATOR)../..$(PATH_SEPARATOR) \
                    WEBMOD=`echo webmod|sed '$(transform)'` \
                    GREP=`echo @GREP@|sed '$(transform)'`


ANDREWS_OUTPUT=\
	$(top_builddir)/Tests/Andrews/webmod.topout \
	$(top_builddir)/Tests/Andrews/andcrk.statvar \
	$(top_builddir)/Tests/Andrews/andcrk.out \
	$(top_builddir)/Tests/Andrews/input/phreeqmms.pqi \
	$(top_builddir)/Tests/Andrews/webmod.chem.out

CLEANFILES=$(ANDREWS_OUTPUT)
