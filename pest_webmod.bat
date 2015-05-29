rm -f andcrk.params
rm -f andcrk.dat.pqi
rm -f andcrk.statvar
par2par par2par_andcrk.dat
/home/rmwebb/programs/webmod-trunk/webmod_1.0.exe -C/home/rmwebb/andcrk_pest3/andcrk.control > webmod.log
/home/dlpark/WEB_PEST/tsproc.exe < /home/rmwebb/andcrk_pest3/tsproc.in
