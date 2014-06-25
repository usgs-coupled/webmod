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
	make --directory="webmod/" --file=Makefile

# Builds project 'mmf_c'...
.PHONY: mmf_c
mmf_c: 
	make --directory="mmf_c/" --file=mmf_c.makefile

# Builds project 'IPhreeqc'...
###.PHONY: IPhreeqc
###IPhreeqc: 
###	make --directory="IPhreeqcMMS/IPhreeqc/" --file=IPhreeqc.makefile

# Builds project 'IPhreeqcMMS'...
.PHONY: IPhreeqcMMS
IPhreeqcMMS: 
	make --directory="IPhreeqcMMS/" --file=Makefile

# Cleans all projects...
.PHONY: clean
clean:
	make --directory="mmf_c/" --file=mmf_c.makefile clean
	make --directory="IPhreeqcMMS/" --file=Makefile clean

