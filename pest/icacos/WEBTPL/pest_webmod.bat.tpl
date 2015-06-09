@DEL@ andcrk.params
@DEL@ andcrk.dat.pqi
@DEL@ andcrk.statvar
@PEST_BIN_DIR@par2par par2par_andcrk.dat
@PROJECT_DIR@webmod_1.0.exe -C@PROJECT_DIR@andcrk.control > webmod.log
@TSPROC_BIN_DIR@tsproc.exe < @PROJECT_DIR@tsproc.in
