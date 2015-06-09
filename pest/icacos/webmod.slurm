#!/bin/sh
#SBATCH -p UV
#SBATCH -A nrp
#SBATCH -n 46
##SBATCH -w compute87
##SBATCH -x compute[188-190,90,148,156,157,159,170,173,174,179,184,194,195]
##SBATCH -x compute[170]
#SBATCH --output=andcrk.pest.out
#SBATCH --mem-per-cpu=2000M
#SBATCH -t72:00:00
##SBATCH --ntasks-per-node=1
##SBATCH --threads-per-core=1
##SBATCH --hint=compute_bound
##SBATCH --hint nomultithread
##SBATCH --exclusive

module load tools/beopest-12.2-gnu
##module load mpi/openmpi-1.5.4
# This should work ...
module load mpi/openmpi-1.6.5-gcc

./clean

TEMP_DIRS=/home/rmwebb/WEB_PEST3
PROJECT_DIR=/home/rmwebb/andcrk_pest3/
WEBMOD_BIN_DIR=/home/rmwebb/programs/webmod-trunk/
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

sed -i "2s/CONTEXT pest_prep/CONTEXT model_run/" ./andcrk_tsproc.dat

# Line 6 of the pst file RLAMBDA1 .....	
sed -i "6s/.*/10.0  -3.0    0.3    0.03     -10  999  LAMFORGIVE/" ./andcrk.pst
# Line 7 of the pst file RELPARMAX FACPARMAX FACORIG.....	
sed -i "7s/.*/0.2   2.0   1.0e-3/" ./andcrk.pst
# Line 9 of the pst file. NOPTMAX Max # of optimizations. Default is 30, set to 0 for single run with phi contributions, 1 for sensitivities, or a small number to test PEST loops.
sed -i "9s/.*/30   .005  4   4  .005   4/" ./andcrk.pst
# Line 13. SVD block MAXSING EIGTHRESH. Replace MAXSING the maximum number of adjustable variables (number of singlular valuess at which truncation occurs)
sed -i "13s/.*/45 5e-7/" ./andcrk.pst
# Line 14 EIGWRITE. 0 if not using SVD output file
sed -i "14s/.*/1/" ./andcrk.pst

# Use this sed copy to t form if run on machine that throws error with -i (in place replace to same file)
# For trial runs, this only does two opimization iterations
#sed "s/30   .005  4   4  .005   4/2   .005  4   4  .005   4/" ./andcrk.pst > t
#mv t ./andcrk.pst

for i in {1..46}
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


