@DEL@ webmod.params
@DEL@ webmod.pqi
@DEL@ webmod.statvar
@PEST_BIN_DIR@par2par par2par.dat
@PROJECT_DIR@webmod.exe -C@PROJECT_DIR@webmod.control > webmod.log
@TSPROC_BIN_DIR@tsproc.exe < @PROJECT_DIR@tsproc.in
