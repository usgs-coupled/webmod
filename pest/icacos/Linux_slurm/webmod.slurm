#!/bin/sh

#SBATCH -p UV
#SBATCH -A nrp
#SBATCH -n 18
#SBATCH --output=ica.pest.out
#SBATCH --mem-per-cpu=2000M
#SBATCH -t72:00:00
#SBATCH --exclusive

module load tools/beopest-12.2-gnu
module load intel/psxe-2015
module load mpi/openmpi-1.6.5-gcc

#./clean

LOCAL_HOME=`pwd`

WEB_DAT=${LOCAL_HOME}/../WEBDATA
WEB_TPL=${LOCAL_HOME}/../WEBTPL
WEB_SRC=./webfiles
PROJECT_DIR=${LOCAL_HOME}/webmod
PEST_BIN_DIR=
TSPROC_BIN_DIR=${LOCAL_HOME}/../../Linux_bin
TEMP_DIRS=${PROJECT_DIR}

# setup WEB_SRC
rm -rf ${WEB_SRC}
mkdir ${WEB_SRC}
sed "s#@PROJECT_DIR@#${PROJECT_DIR}/#" ${WEB_TPL}/ica_tsproc.tpl > ${WEB_SRC}/ica_tsproc.dat
sed "s#@PROJECT_DIR@#${PROJECT_DIR}/#g" ${WEB_TPL}/pest_webmod.bat.tpl > ${WEB_SRC}/pest_webmod.bata
sed "s#@PEST_BIN_DIR@#${PEST_BIN_DIR}#" ${WEB_SRC}/pest_webmod.bata > ${WEB_SRC}/pest_webmod.batb
sed "s#@TSPROC_BIN_DIR@#${TSPROC_BIN_DIR}/#" ${WEB_SRC}/pest_webmod.batb > ${WEB_SRC}/pest_webmod.bat
rm ${WEB_SRC}/pest_webmod.bata ${WEB_SRC}/pest_webmod.batb
sed -i "s#@DEL@#rm -f#" ${WEB_SRC}/pest_webmod.bat

chmod 755 ${WEB_SRC}/pest_webmod.bat
cp ${WEB_DAT}/* ${WEB_SRC}
cp ${TSPROC_BIN_DIR}/webmod_1.0.exe ${WEB_SRC}

# setup working directory PROJECT_DIR
rm -rf ${PROJECT_DIR}
mkdir ${PROJECT_DIR}
cp ${WEB_SRC}/ica.statvar.full ${PROJECT_DIR}/ica.statvar
cp ${WEB_SRC}/* ${PROJECT_DIR}

# run tsproc
cd ${PROJECT_DIR}
${TSPROC_BIN_DIR}/tsproc.exe < ${PROJECT_DIR}/tsproc.in

# Sed CONTEXT ica_tsproc.dat
sed -i "2s/CONTEXT pest_prep/CONTEXT model_run/" ${PROJECT_DIR}/ica_tsproc.dat

# Sed ica.pst
# Line 6 of the pst file RLAMBDA1 .....	
sed -i "6s/.*/10.0  -3.0    0.3    0.03     -10  999  LAMFORGIVE/" ${PROJECT_DIR}/ica.pst
# Line 7 of the pst file RELPARMAX FACPARMAX FACORIG.....	
sed -i "7s/.*/0.2   2.0   1.0e-3/" ${PROJECT_DIR}/ica.pst
# Line 9 of the pst file. NOPTMAX Max # of optimizations. Default is 30, set to 0 for single run with phi contributions, 1 for sensitivities, or a small number to test PEST loops.
sed -i "9s/.*/30   .005  4   4  .005   4/" ${PROJECT_DIR}/ica.pst
# Line 13. SVD block MAXSING EIGTHRESH. Replace MAXSING the maximum number of adjustable variables (number of singlular valuess at which truncation occurs)
sed -i "13s/.*/17 5e-7/" ${PROJECT_DIR}/ica.pst
# Line 14 EIGWRITE. 0 if not using SVD output file
sed -i "14s/.*/1/" ${PROJECT_DIR}/ica.pst

# Make temp directories
for i in {1..18}
do
    rm -rf ica$i
    mkdir ica$i
    cp \
	pqi_ica.tpl \
	params_ica.tpl \
	ica_tsproc.dat \
	ica.ins \
	ica.dat \
	ica.dat.chemdat \
	ica.statvar \
        phreeq_lut \
        phreeqc_web_lite.dat \
	ica$i
done

time mpirun -np 18 ppest ${PROJECT_DIR}/ica /M ${TEMP_DIRS}/ica
