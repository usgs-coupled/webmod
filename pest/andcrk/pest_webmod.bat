rm -f andcrk.params
rm -f andcrk.dat.pqi
rm -f andcrk.statvar
par2par par2par_andcrk.dat
/home/dlpark/WEB_PEST/webmod_1.0.exe -C/home/dlpark/WEB_PEST/andcrk.control > webmod.log
/home/dlpark/WEB_PEST/tsproc.exe < /home/dlpark/WEB_PEST/tsproc.in
