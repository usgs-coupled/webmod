#!/bin/sh
#SBATCH -p Prod
##SBATCH -p Debug
#SBATCH -A Water
#SBATCH -n 20
##SBATCH -w compute87
#SBATCH --output=andcrk.pest.out
#SBATCH --mem-per-cpu=1000M
#SBATCH -t14:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --threads-per-core=1
#SBATCH --hint compute_bound
#SBATCH --hint nomultithread
#SBATCH --exclusive

module load tools/beopest-12.2-gnu
module load mpi/openmpi-1.5.4

./clean

TEMP_DIRS=/home/dlpark/WEB_PEST
PROJECT_DIR=/home/dlpark/andcrk_pest/
WEBMOD_BIN_DIR=/home/dlpark/WEB_PEST/
TSPROC_BIN_DIR=/home/dlpark/WEB_PEST/
PEST_BIN_DIR=

sed "s%@PROJECT_DIR@%${PROJECT_DIR}%" andcrk.control.tpl > andcrk.control
sed "s%@PROJECT_DIR@%${PROJECT_DIR}%" andcrk_tsproc.tpl > andcrk_tsproc.dat
sed "s%@PROJECT_DIR@%${PROJECT_DIR}%" pest_webmod.bat.tpl | \
    sed "s%@WEBMOD_BIN_DIR@%${WEBMOD_BIN_DIR}%" | \
    sed "s%@TSPROC_BIN_DIR@%${TSPROC_BIN_DIR}%" | \
    sed "s%@PEST_BIN_DIR@%${PEST_BIN_DIR}%" > pest_webmod.bat
sed -e "s%@PROJECT_DIR@%${PROJECT_DIR}%" tsproc.in.tpl > tsproc.in
sed -e "s%@PROJECT_DIR@%${PROJECT_DIR}%" par2par_andcrk.tpl.tpl > par2par_andcrk.tpl
cp andcrk.statvar.full andcrk.statvar

${TSPROC_BIN_DIR}tsproc.exe < ./tsproc.in

sed "s/CONTEXT pest_prep/CONTEXT model_run/" ./andcrk_tsproc.dat > t
mv t ./andcrk_tsproc.dat


sed "s/* parameter groups/* singular value decomposition\n1\n4 5e-7\n0\n* parameter groups/" ./andcrk.pst > t
mv t ./andcrk.pst
#sed "s/30   .005  4   4  .005   4/2   .005  4   4  .005   4/" ./andcrk.pst > t
#mv t ./andcrk.pst

for i in {1..20}
do
    rm -rf ${TEMP_DIRS}/andcrk$i
    mkdir ${TEMP_DIRS}/andcrk$i
#    cp * ${TEMP_DIRS}/andcrk$i
    cp \
andcrk.ins \
andcrk.dat \
andcrk.dat.chemdat \
andcrk.statvar \
phreeq_lut \
	${TEMP_DIRS}/andcrk$i
done

time mpirun ${PEST_BIN_DIR}ppest ${PROJECT_DIR}andcrk /M ${TEMP_DIRS}/andcrk


