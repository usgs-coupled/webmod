rm -f andcrk.params
rm -f andcrk.dat.pqi
rm -f andcrk.statvar
@PEST_BIN_DIR@par2par par2par_andcrk.dat
@WEBMOD_BIN_DIR@webmod_1.0.exe -C@PROJECT_DIR@andcrk.control > webmod.log
@TSPROC_BIN_DIR@tsproc.exe < @PROJECT_DIR@tsproc.in
