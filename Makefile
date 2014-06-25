IF_HOME:=/opt/intel/composer_xe_2013.1.117/compiler/lib
IFORT_LIB_64:=-Xlinker -Bstatic -L$(IF_HOME)/intel64/ -lifport -lifcore -limf -lsvml -lm -lipgo -lirc -lirc_s -Xlinker -Bdynamic -lpthread -ldl


LD_FLAGS:=${IFORT_LIB_64}
LINKER:=g++
OBJECT_FILES:=webmod/gccRelease/src/*.o
LIBS:=gccRelease/libmmf_c.a IPhreeqcMMS/lib/libphreeqcmms.a

.PHONY: webm
webm: mmf_c IPhreeqcMMS webmod
	${LINKER} -o $@ ${OBJECT_FILES} ${LIBS} ${LD_FLAGS}

.PHONY: webmod
webmod:
	$(MAKE) -C webmod --file=Makefile

# Builds project 'mmf_c'...
.PHONY: mmf_c
mmf_c: 
	$(MAKE) -C mmf_c --file=mmf_c.makefile

# Builds project 'IPhreeqcMMS'...
.PHONY: IPhreeqcMMS
IPhreeqcMMS: 
	$(MAKE) -C IPhreeqcMMS --file=Makefile

# Cleans all projects...
.PHONY: clean
clean:
	$(MAKE) -C mmf_c --file=mmf_c.makefile clean
	$(MAKE) -C IPhreeqcMMS clean
	$(MAKE) -C webmod --file=Makefile clean

