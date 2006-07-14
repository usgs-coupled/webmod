
FFLAGS     = -Wall -I../include
LDFLAGS    = -L../lib

%.o: %.f
	$(FC) $(FFLAGS) $(TARGET_ARCH) -g -c -o $@ $<

all: kinn

clean:
	$(RM) kinn.o kinn.exe kinn phreeqc.log phreeqc.out phreeqc.err selected.out

kinn: kinn.o ../lib/libiphreeqc.a
	$(CXX) $(LDFLAGS) -o kinn kinn.o -liphreeqc -lfrtbegin -lg2c

kinn.o: kinn.f
