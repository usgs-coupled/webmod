include makelist

CXX        = g++
# CPPFLAGS   = -O3 -Wall -DNDEBUG -DSWIG_SHARED_OBJ
CPPFLAGS   = -g -Wall -DSWIG_SHARED_OBJ
TARGET     = lib/libphreeqcmms.a

TARGET_ARCH = -IIPhreeqc/include

VPATH=src:IPhreeqc/src/phreeqcpp/phreeqc:IPhreeqc/src

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
		dw.o \
		integrate.o \
		input.o \
		inverse.o \
		isotopes.o \
		kinetics.o \
		mainsubs.o \
		model.o \
		nvector.o \
		nvector_serial.o \
		p2clib.o \
		parse.o \
		phqalloc.o \
		pitzer.o \
		pitzer_structures.o \
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
		Phreeqc.o \
		SelectedOutput.o \
		TallyF.o \
		Var.o \
		fwrap.o \
		global.o \
		module_files.o \
		module_output.o


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


# POBJS
advection.o: IPhreeqc/src/phreeqcpp/phreeqc/advection.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
basic.o: IPhreeqc/src/phreeqcpp/phreeqc/basic.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h \
  IPhreeqc/src/phreeqcpp/phreeqc/p2c.h
basicsubs.o: IPhreeqc/src/phreeqcpp/phreeqc/basicsubs.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
cl1.o: IPhreeqc/src/phreeqcpp/phreeqc/cl1.c \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h
cvdense.o: IPhreeqc/src/phreeqcpp/phreeqc/cvdense.c \
  IPhreeqc/src/phreeqcpp/phreeqc/cvdense.h \
  IPhreeqc/src/phreeqcpp/phreeqc/cvode.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/nvector.h \
  IPhreeqc/src/phreeqcpp/phreeqc/dense.h \
  IPhreeqc/src/phreeqcpp/phreeqc/smalldense.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialsmath.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h
cvode.o: IPhreeqc/src/phreeqcpp/phreeqc/cvode.c \
  IPhreeqc/src/phreeqcpp/phreeqc/cvode.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/nvector.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialsmath.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/kinetics.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h
dense.o: IPhreeqc/src/phreeqcpp/phreeqc/dense.c \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialsmath.h \
  IPhreeqc/src/phreeqcpp/phreeqc/dense.h \
  IPhreeqc/src/phreeqcpp/phreeqc/smalldense.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h
dw.o: IPhreeqc/src/phreeqcpp/phreeqc/dw.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/pitzer.h
integrate.o: IPhreeqc/src/phreeqcpp/phreeqc/integrate.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
input.o: IPhreeqc/src/phreeqcpp/phreeqc/input.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/input.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h
inverse.o: IPhreeqc/src/phreeqcpp/phreeqc/inverse.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
isotopes.o: IPhreeqc/src/phreeqcpp/phreeqc/isotopes.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
kinetics.o: IPhreeqc/src/phreeqcpp/phreeqc/kinetics.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqcpp/phreeqc/cvode.h \
  IPhreeqc/src/phreeqcpp/phreeqc/nvector.h \
  IPhreeqc/src/phreeqcpp/phreeqc/cvdense.h \
  IPhreeqc/src/phreeqcpp/phreeqc/dense.h \
  IPhreeqc/src/phreeqcpp/phreeqc/smalldense.h \
  IPhreeqc/src/phreeqcpp/phreeqc/nvector_serial.h \
  IPhreeqc/src/phreeqcpp/phreeqc/kinetics.h
mainsubs.o: IPhreeqc/src/phreeqcpp/phreeqc/mainsubs.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h \
  IPhreeqc/src/phreeqcpp/phreeqc/input.h
model.o: IPhreeqc/src/phreeqcpp/phreeqc/model.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
nvector.o: IPhreeqc/src/phreeqcpp/phreeqc/nvector.c \
  IPhreeqc/src/phreeqcpp/phreeqc/nvector.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h
nvector_serial.o: IPhreeqc/src/phreeqcpp/phreeqc/nvector_serial.c \
  IPhreeqc/src/phreeqcpp/phreeqc/nvector_serial.h \
  IPhreeqc/src/phreeqcpp/phreeqc/nvector.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialsmath.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h
p2clib.o: IPhreeqc/src/phreeqcpp/phreeqc/p2clib.c \
  IPhreeqc/src/phreeqcpp/phreeqc/p2c.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h
parse.o: IPhreeqc/src/phreeqcpp/phreeqc/parse.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
phqalloc.o: IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h
pitzer.o: IPhreeqc/src/phreeqcpp/phreeqc/pitzer.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h \
  IPhreeqc/src/phreeqcpp/phreeqc/pitzer.h
pitzer_structures.o: IPhreeqc/src/phreeqcpp/phreeqc/pitzer_structures.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h \
  IPhreeqc/src/phreeqcpp/phreeqc/pitzer.h
prep.o: IPhreeqc/src/phreeqcpp/phreeqc/prep.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
print.o: IPhreeqc/src/phreeqcpp/phreeqc/print.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h \
  IPhreeqc/src/phreeqcpp/phreeqc/pitzer.h
read.o: IPhreeqc/src/phreeqcpp/phreeqc/read.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
readtr.o: IPhreeqc/src/phreeqcpp/phreeqc/readtr.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
smalldense.o: IPhreeqc/src/phreeqcpp/phreeqc/smalldense.c \
  IPhreeqc/src/phreeqcpp/phreeqc/smalldense.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialsmath.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h
spread.o: IPhreeqc/src/phreeqcpp/phreeqc/spread.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
step.o: IPhreeqc/src/phreeqcpp/phreeqc/step.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
structures.o: IPhreeqc/src/phreeqcpp/phreeqc/structures.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
sundialsmath.o: IPhreeqc/src/phreeqcpp/phreeqc/sundialsmath.c \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialsmath.h \
  IPhreeqc/src/phreeqcpp/phreeqc/sundialstypes.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h
tally.o: IPhreeqc/src/phreeqcpp/phreeqc/tally.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
tidy.o: IPhreeqc/src/phreeqcpp/phreeqc/tidy.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
transport.o: IPhreeqc/src/phreeqcpp/phreeqc/transport.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
utilities.o: IPhreeqc/src/phreeqcpp/phreeqc/utilities.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h

# SOBJS
IPhreeqc.o: IPhreeqc/src/IPhreeqc.cxx IPhreeqc/src/phreeqcns.hxx \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/input.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h IPhreeqc/src/Phreeqc.hxx \
  IPhreeqc/src/ErrorReporter.hxx IPhreeqc/src/SelectedOutput.hxx \
  IPhreeqc/src/CVar.hxx IPhreeqc/src/Debug.h \
  IPhreeqc/src/../include/Var.h IPhreeqc/src/../include/IPhreeqc.h \
  IPhreeqc/src/../include/Var.h IPhreeqc/src/module_files.h
IPhreeqcF.o: IPhreeqc/src/IPhreeqcF.F
Phreeqc.o: IPhreeqc/src/Phreeqc.cxx
SelectedOutput.o: IPhreeqc/src/SelectedOutput.cxx \
  IPhreeqc/src/SelectedOutput.hxx IPhreeqc/src/CVar.hxx \
  IPhreeqc/src/Debug.h IPhreeqc/src/../include/Var.h \
  IPhreeqc/src/phreeqcns.hxx IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/input.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h
TallyF.o: IPhreeqc/src/TallyF.F
Var.o: IPhreeqc/src/Var.c IPhreeqc/src/../include/Var.h
fwrap.o: IPhreeqc/src/fwrap.c IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/../include/IPhreeqc.h IPhreeqc/src/../include/Var.h
module_files.o: IPhreeqc/src/module_files.c IPhreeqc/src/module_files.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phreeqc_files.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h \
  IPhreeqc/src/phreeqcpp/phreeqc/input.h
module_output.o: IPhreeqc/src/module_output.c IPhreeqc/src/module_files.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.c \
  IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  IPhreeqc/src/phreeqcpp/phreeqc/output.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phrqproto.h \
  IPhreeqc/src/phreeqcpp/phreeqc/phqalloc.h

# MOBJS
phr_cmix.o: src/phr_cmix.c src/../IPhreeqc/src/phreeqcpp/phreeqc/global.h \
  src/../IPhreeqc/src/phreeqcpp/phreeqc/phrqtype.h \
  src/../IPhreeqc/include/IPhreeqc.h src/../IPhreeqc/include/Var.h
phr_mix.o: src/phr_mix.F
phr_multicopy.o: src/phr_multicopy.f
phr_precip.o: src/phr_precip.f