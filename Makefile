include makelist

CXX        = g++
# CPPFLAGS   = -O3 -Wall -DNDEBUG -DSWIG_SHARED_OBJ
CPPFLAGS   = -g -Wall -DSWIG_SHARED_OBJ
TARGET     = lib/libphreeqcmms.a

VPATH=src:IPhreeqc/src/phreeqc:IPhreeqc/src

%.o: %.f
	$(FC) $(FFLAGS) $(TARGET_ARCH) -c -o $@ $<

%.o: %.F
	$(FC) $(FFLAGS) $(TARGET_ARCH) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c -o $@ $<

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c -o $@ $<

PDIR = 	IPhreeqc/src/phreeqc
POBJS =	\
		advection.o \
		basic.o \
		basicsubs.o \
		cl1.o \
		cvdense.o \
		cvode.o \
		dense.o \
		integrate.o \
		inverse.o \
		isotopes.o \
		kinetics.o \
		mainsubs.o \
		message.o \
		model.o \
		nvector.o \
		nvector_serial.o \
		p2clib.o \
		parse.o \
		phqalloc.o \
		prep.o \
		print.o \
		read.o \
		readtr.o \
		smalldense.o \
		spread.o \
		step.o \
		structures.o \
		sundialsmath.o \
		tally.o \
		tidy.o \
		transport.o \
		utilities.o


SOBJS =	\
		IPhreeqc.o \
		IPhreeqcF.o \
		Output.o \
		Overrides.o \
		Phreeqc.o \
		PhreeqcParser.o \
		SelectedOutput.o \
		TallyF.o \
		Var.o \
		fwrap.o \
		global.o


MOBJS =	\
		phr_cmix.o \
		phr_mix.o \
		phr_multicopy.o \
		phr_precip.o

all: $(TARGET)


$(TARGET): $(POBJS) $(SOBJS) $(MOBJS) lib
	$(AR)  $(TARGET) $(POBJS) $(SOBJS) $(MOBJS)
	$(RANLIB) $(TARGET)	

lib:
	mkdir -p lib

clean:
	$(RM) $(POBJS) $(SOBJS) $(MOBJS)



# iphreeqcmms (MOBJS)
phr_cmix.o: src/phr_cmix.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/include/IPhreeqc.h IPhreeqc/include/Var.h
phr_mix.o: src/phr_mix.F IPhreeqc/include/IPhreeqc.f.inc
phr_multicopy.o: src/phr_multicopy.f  IPhreeqc/include/IPhreeqc.f.inc
phr_precip.o: src/phr_precip.f  IPhreeqc/include/IPhreeqc.f.inc

# iphreeqc (SOBJS)
# C++
IPhreeqc.o: IPhreeqc/src/IPhreeqc.cxx IPhreeqc/src/phreeqcns.hxx \
  IPhreeqc/src/phreeqc/global.h IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h IPhreeqc/src/Phreeqc.hxx \
  IPhreeqc/include/IPhreeqc.h IPhreeqc/include/Var.h \
  IPhreeqc/src/ErrorReporter.hxx IPhreeqc/src/SelectedOutput.hxx \
  IPhreeqc/src/CVar.hxx IPhreeqc/src/Debug.h
Output.o: IPhreeqc/src/Output.cxx IPhreeqc/src/Output.hxx \
  IPhreeqc/src/Debug.h IPhreeqc/src/phreeqcns.hxx \
  IPhreeqc/src/phreeqc/global.h IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h IPhreeqc/src/ErrorReporter.hxx
Overrides.o: IPhreeqc/src/Overrides.cxx IPhreeqc/src/Debug.h \
  IPhreeqc/src/Phreeqc.hxx IPhreeqc/include/IPhreeqc.h \
  IPhreeqc/include/Var.h IPhreeqc/src/PhreeqcParser.hxx
Phreeqc.o: IPhreeqc/src/Phreeqc.cxx IPhreeqc/src/Debug.h \
  IPhreeqc/src/Phreeqc.hxx IPhreeqc/include/IPhreeqc.h \
  IPhreeqc/include/Var.h IPhreeqc/src/phreeqcns.hxx \
  IPhreeqc/src/phreeqc/global.h IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h IPhreeqc/src/PhreeqcParser.hxx \
  IPhreeqc/src/Output.hxx IPhreeqc/src/ErrorReporter.hxx
PhreeqcParser.o: IPhreeqc/src/PhreeqcParser.cxx \
  IPhreeqc/src/PhreeqcParser.hxx IPhreeqc/src/phreeqcns.hxx \
  IPhreeqc/src/phreeqc/global.h IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h IPhreeqc/src/Debug.h
SelectedOutput.o: IPhreeqc/src/SelectedOutput.cxx \
  IPhreeqc/src/SelectedOutput.hxx IPhreeqc/src/CVar.hxx \
  IPhreeqc/src/Debug.h IPhreeqc/include/Var.h IPhreeqc/src/phreeqcns.hxx \
  IPhreeqc/src/phreeqc/global.h IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h
# Fortran
IPhreeqcF.o: IPhreeqc/src/IPhreeqcF.F
TallyF.o: IPhreeqc/src/TallyF.F
# C
Var.o: IPhreeqc/src/Var.c IPhreeqc/include/Var.h
fwrap.o: IPhreeqc/src/fwrap.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/include/IPhreeqc.h IPhreeqc/include/Var.h
global.o: IPhreeqc/src/global.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h

Var.o: IPhreeqc/src/Var.c IPhreeqc/include/Var.h
fwrap.o: IPhreeqc/src/fwrap.c IPhreeqc/include/IPhreeqc.h \
  IPhreeqc/include/Var.h
global.o: IPhreeqc/src/global.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h

# phreeqc (POBJS)
advection.o: IPhreeqc/src/phreeqc/advection.c \
  IPhreeqc/src/phreeqc/global.h IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h
basic.o: IPhreeqc/src/phreeqc/basic.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h \
  IPhreeqc/src/phreeqc/p2c.h
basicsubs.o: IPhreeqc/src/phreeqc/basicsubs.c \
  IPhreeqc/src/phreeqc/global.h IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h
cl1.o: IPhreeqc/src/phreeqc/cl1.c IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h
cvdense.o: IPhreeqc/src/phreeqc/cvdense.c IPhreeqc/src/phreeqc/cvdense.h \
  IPhreeqc/src/phreeqc/cvode.h IPhreeqc/src/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqc/nvector.h IPhreeqc/src/phreeqc/dense.h \
  IPhreeqc/src/phreeqc/smalldense.h IPhreeqc/src/phreeqc/sundialsmath.h \
  IPhreeqc/src/phreeqc/message.h
cvode.o: IPhreeqc/src/phreeqc/cvode.c IPhreeqc/src/phreeqc/cvode.h \
  IPhreeqc/src/phreeqc/sundialstypes.h IPhreeqc/src/phreeqc/nvector.h \
  IPhreeqc/src/phreeqc/sundialsmath.h IPhreeqc/src/phreeqc/message.h \
  IPhreeqc/src/phreeqc/kinetics.h
dense.o: IPhreeqc/src/phreeqc/dense.c \
  IPhreeqc/src/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqc/sundialsmath.h IPhreeqc/src/phreeqc/dense.h \
  IPhreeqc/src/phreeqc/smalldense.h IPhreeqc/src/phreeqc/message.h
integrate.o: IPhreeqc/src/phreeqc/integrate.c \
  IPhreeqc/src/phreeqc/global.h IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h
inverse.o: IPhreeqc/src/phreeqc/inverse.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
isotopes.o: IPhreeqc/src/phreeqc/isotopes.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
kinetics.o: IPhreeqc/src/phreeqc/kinetics.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h \
  IPhreeqc/src/phreeqc/sundialstypes.h IPhreeqc/src/phreeqc/cvode.h \
  IPhreeqc/src/phreeqc/nvector.h IPhreeqc/src/phreeqc/cvdense.h \
  IPhreeqc/src/phreeqc/dense.h IPhreeqc/src/phreeqc/smalldense.h \
  IPhreeqc/src/phreeqc/nvector_serial.h IPhreeqc/src/phreeqc/kinetics.h
main.o: IPhreeqc/src/phreeqc/main.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
mainsubs.o: IPhreeqc/src/phreeqc/mainsubs.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
message.o: IPhreeqc/src/phreeqc/message.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
model.o: IPhreeqc/src/phreeqc/model.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
nvector.o: IPhreeqc/src/phreeqc/nvector.c IPhreeqc/src/phreeqc/nvector.h \
  IPhreeqc/src/phreeqc/sundialstypes.h IPhreeqc/src/phreeqc/message.h
nvector_serial.o: IPhreeqc/src/phreeqc/nvector_serial.c \
  IPhreeqc/src/phreeqc/nvector_serial.h IPhreeqc/src/phreeqc/nvector.h \
  IPhreeqc/src/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqc/sundialsmath.h IPhreeqc/src/phreeqc/message.h
p2clib.o: IPhreeqc/src/phreeqc/p2clib.c IPhreeqc/src/phreeqc/p2c.h \
  IPhreeqc/src/phreeqc/message.h
parse.o: IPhreeqc/src/phreeqc/parse.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
phqalloc.o: IPhreeqc/src/phreeqc/phqalloc.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/message.h
prep.o: IPhreeqc/src/phreeqc/prep.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
print.o: IPhreeqc/src/phreeqc/print.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
read.o: IPhreeqc/src/phreeqc/read.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
readtr.o: IPhreeqc/src/phreeqc/readtr.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
smalldense.o: IPhreeqc/src/phreeqc/smalldense.c \
  IPhreeqc/src/phreeqc/smalldense.h IPhreeqc/src/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqc/sundialsmath.h IPhreeqc/src/phreeqc/message.h
spread.o: IPhreeqc/src/phreeqc/spread.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
step.o: IPhreeqc/src/phreeqc/step.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
structures.o: IPhreeqc/src/phreeqc/structures.c \
  IPhreeqc/src/phreeqc/global.h IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h
sundialsmath.o: IPhreeqc/src/phreeqc/sundialsmath.c \
  IPhreeqc/src/phreeqc/sundialsmath.h \
  IPhreeqc/src/phreeqc/sundialstypes.h IPhreeqc/src/phreeqc/message.h
tally.o: IPhreeqc/src/phreeqc/tally.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
tidy.o: IPhreeqc/src/phreeqc/tidy.c IPhreeqc/src/phreeqc/global.h \
  IPhreeqc/src/phreeqc/phqalloc.h IPhreeqc/src/phreeqc/message.h
transport.o: IPhreeqc/src/phreeqc/transport.c \
  IPhreeqc/src/phreeqc/global.h IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h
utilities.o: IPhreeqc/src/phreeqc/utilities.c \
  IPhreeqc/src/phreeqc/global.h IPhreeqc/src/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqc/message.h