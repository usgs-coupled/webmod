# Builds all the projects in the solution...
.PHONY: all_projects
all_projects: mmf_c IPhreeqc IPhreeqcMMS 

# Builds project 'mmf_c'...
.PHONY: mmf_c
mmf_c: 
	make --directory="mmf_c/" --file=mmf_c.makefile

# Builds project 'IPhreeqc'...
.PHONY: IPhreeqc
IPhreeqc: 
	make --directory="IPhreeqcMMS/IPhreeqc/" --file=IPhreeqc.makefile

# Builds project 'IPhreeqcMMS'...
.PHONY: IPhreeqcMMS
IPhreeqcMMS: 
	make --directory="IPhreeqcMMS/" --file=IPhreeqcMMS.makefile

# Cleans all projects...
.PHONY: clean
clean:
	make --directory="mmf_c/" --file=mmf_c.makefile clean
	make --directory="IPhreeqcMMS/IPhreeqc/" --file=IPhreeqc.makefile clean
	make --directory="IPhreeqcMMS/" --file=IPhreeqcMMS.makefile clean

