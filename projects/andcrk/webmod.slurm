#!/bin/sh

PROCESSES=5
#SBATCH -p prod
#SBATCH -A nrp
#SBATCH --output=webmod.pest.out
#SBATCH --mem-per-cpu=1000M
#SBATCH -t00:10:00
#SBATCH --exclusive

module load tools/beopest-12.2-gnu
module load intel/psxe-2015
module load mpi/openmpi-1.6.5-gcc

LOCAL_HOME=`pwd`
INPUT_DIR=${LOCAL_HOME}/input
PEST_FILES_DIR=${LOCAL_HOME}/pest_files
CONTROL_DIR=${LOCAL_HOME}/control
PROJECT_DIR=${LOCAL_HOME}/pest_run_dir
BIN_DIR=${LOCAL_HOME}/../../bin
PEST_BIN_DIR=
TEMP_DIR=psttemp

# setup working directory PROJECT_DIR=pest_run_dir
rm -rf ${PROJECT_DIR}
mkdir ${PROJECT_DIR}

# Sed to make tsproc.dat 
sed "s#@PROJECT_DIR@#${PROJECT_DIR}/#g"    ${PEST_FILES_DIR}/tsproc.dat.tpl > ${PROJECT_DIR}/tsproc.dat

#Sed to make pest_webmod.bat in PROJECT_DIR
sed "s#@PROJECT_DIR@#${PROJECT_DIR}/#g"    ${PEST_FILES_DIR}/pest_webmod.bat.tpl >  ${PROJECT_DIR}/pest_webmod.bat
sed -i "s#@PEST_BIN_DIR@#${PEST_BIN_DIR}#" ${PROJECT_DIR}/pest_webmod.bat
sed -i "s#@TSPROC_BIN_DIR@#${BIN_DIR}/#"   ${PROJECT_DIR}/pest_webmod.bat
sed -i "s#@DEL@#rm -f#"                    ${PROJECT_DIR}/pest_webmod.bat
sed -i "s#\.exe##g"                        ${PROJECT_DIR}/pest_webmod.bat
chmod 755 ${PROJECT_DIR}/pest_webmod.bat

# Copy files to tsproc directory
cp ${PEST_FILES_DIR}/*.ssf            	   ${PROJECT_DIR}
cp ${PEST_FILES_DIR}/pest_groups*.txt 	   ${PROJECT_DIR}
cp ${PEST_FILES_DIR}/pest_params*.txt 	   ${PROJECT_DIR}
cp ${PEST_FILES_DIR}/tsproc.in         	   ${PROJECT_DIR}
cp ${PEST_FILES_DIR}/webmod.statvar	   ${PROJECT_DIR}
cp ${PEST_FILES_DIR}/par2par.dat.tpl	   ${PROJECT_DIR}

# run tsproc
cd ${PROJECT_DIR}
${BIN_DIR}/tsproc < ${PROJECT_DIR}/tsproc.in

# edit the control file
sed "s#\./input#\.#"      ${CONTROL_DIR}/webmod.control  >  ${PROJECT_DIR}/webmod.control
sed -i "s#\./output#\.#"  ${PROJECT_DIR}/webmod.control 
sed -i "s#\.\./\.\./bin#${BIN_DIR}#"  ${PROJECT_DIR}/webmod.control 

# copy executable
cp ${BIN_DIR}/webmod ${PROJECT_DIR}/webmod

# Sed CONTEXT tsproc.dat
sed -i "2s/CONTEXT pest_prep/CONTEXT model_run/" ${PROJECT_DIR}/tsproc.dat

# Sed andcrk.pst
# Line 6 of the pst file RLAMBDA1 .....	
sed -i "6s/.*/10.0  -3.0    0.3    0.03     -10  999  LAMFORGIVE/" ${PROJECT_DIR}/webmod.pst
# Line 7 of the pst file RELPARMAX FACPARMAX FACORIG.....	
sed -i "7s/.*/0.2   2.0   1.0e-3/" ${PROJECT_DIR}/webmod.pst
# Line 9 of the pst file. NOPTMAX Max # of optimizations. Default is 30, set to 0 for single run with phi contributions, 1 for sensitivities, or a small number to test PEST loops.
sed -i "9s/.*/30   .005  4   4  .005   4/" ${PROJECT_DIR}/webmod.pst
# Line 13. SVD block MAXSING EIGTHRESH. Replace MAXSING the maximum number of adjustable variables (number of singlular valuess at which truncation occurs)
sed -i "13s/.*/19 5e-7/" ${PROJECT_DIR}/webmod.pst
# Line 14 EIGWRITE. 0 if not using SVD output file
sed -i "14s/.*/1/" ${PROJECT_DIR}/webmod.pst

# Make temp directories
i=1
while [ "$i" -le "${PROCESSES}" ]; do
#for i in {1.."${PROCESSES}"}
#do
    cd ${PROJECT_DIR}
    rm -rf ${TEMP_DIR}$i
    mkdir ${TEMP_DIR}$i
    cd ${TEMP_DIR}$i
    cp ${PEST_FILES_DIR}/webmod.params.tpl .
    cp ${PEST_FILES_DIR}/webmod.pqi.tpl    .
    cp ${INPUT_DIR}/webmod.hydro.dat       .
    cp ${INPUT_DIR}/webmod.chem.dat        .
    cp ${INPUT_DIR}/phreeq_lut             .
    cp ${INPUT_DIR}/phreeqc_web_lite.dat   .
    cp ${PROJECT_DIR}/pest_webmod.bat      .
    cp ${PROJECT_DIR}/tsproc.sim.out.ins   .
    cp ${PROJECT_DIR}/tsproc.dat           .
    i=$(($i+1))
done

# /L requires a run in PROJECT_DIR
cd ${PROJECT_DIR}
cp ${PEST_FILES_DIR}/webmod.params.tpl .
cp ${PEST_FILES_DIR}/webmod.pqi.tpl    .
cp ${INPUT_DIR}/webmod.hydro.dat       .
cp ${INPUT_DIR}/webmod.chem.dat        .
cp ${INPUT_DIR}/phreeq_lut             .
cp ${INPUT_DIR}/phreeqc_web_lite.dat   .

mpirun -np ${PROCESSES} ppest ${PROJECT_DIR}/webmod.pst /M /L ${PROJECT_DIR}/psttemp
#mpirun -np 5 ppest ${PROJECT_DIR}/webmod.pst /M 

