@DEL@ ica.params
@DEL@ ica.dat.pqi
@DEL@ ica.statvar
@PEST_BIN_DIR@par2par par2par_ica.dat
@PROJECT_DIR@webmod_1.0.exe -C@PROJECT_DIR@ica.control > webmod.log
@TSPROC_BIN_DIR@tsproc.exe < @PROJECT_DIR@tsproc.in
