del andcrk.params
del andcrk.dat.pqi
del andcrk.statvar
@PEST_BIN_DIR@par2par par2par_andcrk.dat
@PROJECT_DIR@webmod_1.0.exe -C@PROJECT_DIR@andcrk.control > webmod.log
@PEST_BIN_DIR@tsproc.exe < @PROJECT_DIR@tsproc.in
