FORTRAN_SOURCES+=$(shell find webmod.d/src -name "*.f90" -o -name "*.F90" -o -name "*.f" -o -name "*.F" | egrep -v "/Tests/" | egrep -v "module_dummy" | sed 's^\./^^g')

FORTRAN_SOURCES+=$(shell find IPhreeqcMMS/src -name "*.f90" -o -name "*.F90" -o -name "*.f" -o -name "*.F" | egrep -v "/Tests/" | egrep -v "module_dummy" | sed 's^\./^^g')

FORTRAN_SOURCES+=$(shell find IPhreeqcMMS/IPhreeqc/src -name "*.f90" -o -name "*.F90" -o -name "*.f" -o -name "*.F" | egrep -v "/Tests/" | egrep -v "module_dummy" | sed 's^\./^^g')

depend .depend:
	@makedepf90 $(FORTRAN_SOURCES) | sed 's/\.o/\.\$$\(OBJEXT\)/g' | sort > .depend
