CXX        = g++
CPPFLAGS   = -Wall -I./clib -DNDEBUG -DSWIG_SHARED_OBJ

%.o: %.f
	$(FC) $(FFLAGS) $(TARGET_ARCH) -g -c -o $@ $<

all: test2

clean:
	$(RM) test2.o test2.exe phreeqc.log phreeqc.out phreeqc.err error.inp selected.out

test2: test2.o ../lib/libiphreeqc.a
	$(CXX) $(LDFLAGS) -o test2 -L../lib test2.o -liphreeqc -lfrtbegin -lg2c

test2.o: test2.f
