CC                  = gcc
CFLAGS              = -Wall
##CPPFLAGS            = -O3 -DNDEBUG -DSWIG_SHARED_OBJ
CPPFLAGS            = -g -DSWIG_SHARED_OBJ  # debug
TARGET_ARCH         =
FC                  = g77
FFLAGS              = -fno-second-underscore
CXX                 = g++
CXXFLAGS            = -Wall
RANLIB              = ranlib
RM                  = rm -f
AR                  = ar ruv
TARGET              = lib/libiphreeqc.a

VPATH=src:src/phreeqc

%.o: %.f
	$(FC) $(FFLAGS) $(TARGET_ARCH) -c -o $@ $<

%.o: %.F
	$(FC) $(FFLAGS) $(TARGET_ARCH) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c -o $@ $<

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c -o $@ $<

	
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
		Var.o \
		fwrap.o \
		global.o

all: $(TARGET)

$(TARGET): $(POBJS) $(SOBJS)
	$(AR)  $(TARGET) $(POBJS) $(SOBJS)
	$(RANLIB) $(TARGET)	

clean:
	$(RM) $(POBJS) $(SOBJS)


# iphreeqc (SOBJS)
IPhreeqc.o: src/IPhreeqc.cxx src/phreeqcns.hxx src/phreeqc/global.h \
  src/phreeqc/phqalloc.h src/Phreeqc.hxx src/ErrorReporter.hxx \
  src/SelectedOutput.hxx src/CVar.hxx src/Debug.h include/Var.h \
  include/IPhreeqc.h
Output.o: src/Output.cxx src/Output.hxx src/Debug.h src/phreeqcns.hxx \
  src/phreeqc/global.h src/phreeqc/phqalloc.h src/ErrorReporter.hxx
Overrides.o: src/Overrides.cxx src/Debug.h src/Phreeqc.hxx \
  src/PhreeqcParser.hxx
Phreeqc.o: src/Phreeqc.cxx src/Debug.h src/Phreeqc.hxx src/phreeqcns.hxx \
  src/phreeqc/global.h src/phreeqc/phqalloc.h src/PhreeqcParser.hxx \
  src/Output.hxx src/ErrorReporter.hxx
PhreeqcParser.o: src/PhreeqcParser.cxx src/PhreeqcParser.hxx \
  src/phreeqcns.hxx src/phreeqc/global.h src/phreeqc/phqalloc.h \
  src/Debug.h
SelectedOutput.o: src/SelectedOutput.cxx src/SelectedOutput.hxx \
  src/CVar.hxx src/Debug.h include/Var.h src/phreeqcns.hxx \
  src/phreeqc/global.h src/phreeqc/phqalloc.h
Var.o: src/Var.c include/Var.h
fwrap.o: src/fwrap.c include/IPhreeqc.h include/Var.h
global.o: src/global.c src/phreeqc/global.h src/phreeqc/phqalloc.h
IPhreeqcF.o: src/IPhreeqcF.F


# phreeqc (POBJS)
advection.o: src/phreeqc/advection.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h
basic.o: src/phreeqc/basic.c src/phreeqc/global.h src/phreeqc/phqalloc.h \
  src/phreeqc/p2c.h
basicsubs.o: src/phreeqc/basicsubs.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h
cl1.o: src/phreeqc/cl1.c src/phreeqc/phqalloc.h
cvdense.o: src/phreeqc/cvdense.c src/phreeqc/cvdense.h \
  src/phreeqc/cvode.h src/phreeqc/sundialstypes.h src/phreeqc/nvector.h \
  src/phreeqc/dense.h src/phreeqc/smalldense.h src/phreeqc/sundialsmath.h
cvode.o: src/phreeqc/cvode.c src/phreeqc/cvode.h \
  src/phreeqc/sundialstypes.h src/phreeqc/nvector.h \
  src/phreeqc/sundialsmath.h src/phreeqc/kinetics.h
dense.o: src/phreeqc/dense.c src/phreeqc/sundialstypes.h \
  src/phreeqc/sundialsmath.h src/phreeqc/dense.h src/phreeqc/smalldense.h
integrate.o: src/phreeqc/integrate.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h
inverse.o: src/phreeqc/inverse.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h
isotopes.o: src/phreeqc/isotopes.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h
kinetics.o: src/phreeqc/kinetics.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h src/phreeqc/sundialstypes.h src/phreeqc/cvode.h \
  src/phreeqc/nvector.h src/phreeqc/cvdense.h src/phreeqc/dense.h \
  src/phreeqc/smalldense.h src/phreeqc/nvector_serial.h \
  src/phreeqc/kinetics.h
mainsubs.o: src/phreeqc/mainsubs.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h
message.o: src/phreeqc/message.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h src/phreeqc/message.h
model.o: src/phreeqc/model.c src/phreeqc/global.h src/phreeqc/phqalloc.h
nvector.o: src/phreeqc/nvector.c src/phreeqc/nvector.h \
  src/phreeqc/sundialstypes.h
nvector_serial.o: src/phreeqc/nvector_serial.c \
  src/phreeqc/nvector_serial.h src/phreeqc/nvector.h \
  src/phreeqc/sundialstypes.h src/phreeqc/sundialsmath.h
p2clib.o: src/phreeqc/p2clib.c src/phreeqc/p2c.h
parse.o: src/phreeqc/parse.c src/phreeqc/global.h src/phreeqc/phqalloc.h
phqalloc.o: src/phreeqc/phqalloc.c src/phreeqc/global.h
prep.o: src/phreeqc/prep.c src/phreeqc/global.h src/phreeqc/phqalloc.h
print.o: src/phreeqc/print.c src/phreeqc/global.h src/phreeqc/phqalloc.h
read.o: src/phreeqc/read.c src/phreeqc/global.h src/phreeqc/phqalloc.h
readtr.o: src/phreeqc/readtr.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h
smalldense.o: src/phreeqc/smalldense.c src/phreeqc/smalldense.h \
  src/phreeqc/sundialstypes.h src/phreeqc/sundialsmath.h
spread.o: src/phreeqc/spread.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h
step.o: src/phreeqc/step.c src/phreeqc/global.h src/phreeqc/phqalloc.h
structures.o: src/phreeqc/structures.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h
sundialsmath.o: src/phreeqc/sundialsmath.c src/phreeqc/sundialsmath.h \
  src/phreeqc/sundialstypes.h
tally.o: src/phreeqc/tally.c src/phreeqc/global.h src/phreeqc/phqalloc.h \
  src/phreeqc/message.h
tidy.o: src/phreeqc/tidy.c src/phreeqc/global.h src/phreeqc/phqalloc.h
transport.o: src/phreeqc/transport.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h
utilities.o: src/phreeqc/utilities.c src/phreeqc/global.h \
  src/phreeqc/phqalloc.h


