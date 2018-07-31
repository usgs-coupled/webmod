include ../makelist

CXX        = g++
CPPFLAGS   = -Wall -I./clib -DNDEBUG -DSWIG_SHARED_OBJ

%.o: %.f
	$(FC) $(FFLAGS) $(TARGET_ARCH) -g -c -o $@ $<

all: webb

clean:
	$(RM) webb.o webb webb.exe sel1 sel2 phreeqc.log phreeqc.out phreeqc.err

webb: webb.o ../lib/libphreeqcmms.a
	$(CXX) $(LDFLAGS) -o webb -L../lib webb.o -lphreeqcmms -lfrtbegin -lg2c

webb.o: webb.f
