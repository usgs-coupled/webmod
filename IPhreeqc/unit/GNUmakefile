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
TARGET              = unit
LDFLAGS             = -L/usr/local/lib

VPATH=.:../src

%.o: %.f
	$(FC) $(FFLAGS) $(TARGET_ARCH) -c -o $@ $<

%.o: %.F
	$(FC) $(FFLAGS) $(TARGET_ARCH) -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c -o $@ $<

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -c -o $@ $<

OBJS =	\
		SelectedOutput.o \
		Var.o \
		TestCVar.o \
		TestSelectedOutput.o \
		TestVar.o \
		unit.o 

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(LDFLAGS) -o unit $(OBJS) -lcppunit

clean:
	$(RM) $(OBJS)


# OBJS
SelectedOutput.o: ../src/SelectedOutput.cxx ../src/SelectedOutput.hxx \
  ../src/CVar.hxx ../src/Debug.h ../include/Var.h ../src/phreeqcns.hxx \
  ../src/phreeqc/global.h ../src/phreeqc/phqalloc.h
Var.o: ../src/Var.c ../include/Var.h
TestCVar.o: TestCVar.cpp TestCVar.h ../src/CVar.hxx ../src/Debug.h \
  ../include/Var.h
TestSelectedOutput.o: TestSelectedOutput.cpp TestSelectedOutput.h \
  ../src/SelectedOutput.hxx ../src/CVar.hxx ../src/Debug.h \
  ../include/Var.h
TestVar.o: TestVar.cpp TestVar.h ../include/Var.h
unit.o: unit.cpp TestVar.h ../include/Var.h TestCVar.h ../src/CVar.hxx \
  ../src/Debug.h TestSelectedOutput.h ../src/SelectedOutput.hxx


