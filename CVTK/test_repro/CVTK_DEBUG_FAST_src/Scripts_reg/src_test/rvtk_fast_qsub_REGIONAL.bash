#!/bin/bash
#
#======================================================================
# CROCO is a branch of ROMS developped at IRD and INRIA, in France
# The two other branches from UCLA (Shchepetkin et al) 
# and Rutgers University (Arango et al) are under MIT/X style license.
# CROCO specific routines (nesting) are under CeCILL-C license.
# 
# CROCO website : http://www.croco-ocean.org
#======================================================================
#
#---------------------------------------------------------------------
# Script to Run CVTK DEBUG procedure managing parallelization type 
# AND AGRIF nesting type (No nesting, Nesting 1-way, Nesting 2-ways) : 
# VORTEX and REGIONAL
#--------------------------------------------------------------------
#set -x

echo "=============================================="
echo "=> CONFIG "$mytest

echo "Remove *.exe* *.log* "
[ ! -z "$(ls *.exe* 2>/dev/null)" ] && /bin/rm *.exe*
[ ! -z "$(ls *.log* 2>/dev/null)" ] &&/bin/rm *.log*
echo "Remove the CHECKFILE"
[ -f check_file ] && /bin/rm check_file

#echo "Remove AGRIF_FixedGrids.in"
#/bin/rm -f AGRIF_FixedGrids.in 
#echo " " 

#=============================================================================================
#===
source CONFIGURE_GLOBAL
#===
echo " "
echo "=> MPIRUN COMMAND: "$MPIRUN
#
SOURCE_CVTK=${SOURCE_CROCO}/../CVTK/test_repro/CVTK_DEBUG_FAST_src
echo 'Sources CVTK tests: '$SOURCE_CVTK
#
# Get updated files
#
\cp ${SOURCE_CROCO}/cppdefs.h cppdefs.h.SERIAL
sed 's/'undef\ \ \*RVTK_DEBUG'/'define\ RVTK_DEBUG'/' < cppdefs.h.SERIAL > cppdefs.h.SERIAL.tmp
mv cppdefs.h.SERIAL.tmp cppdefs.h.SERIAL
\cp ${SOURCE_CROCO}/param.h param.h.SERIAL

\cp ${SOURCE_CROCO}/cppdefs.h cppdefs.h.OPENMP
sed 's/'undef\ \ \*RVTK_DEBUG'/'define\ RVTK_DEBUG'/' < cppdefs.h.OPENMP > cppdefs.h.OPENMP.tmp
mv cppdefs.h.OPENMP.tmp cppdefs.h.OPENMP
\cp ${SOURCE_CROCO}/param.h param.h.OPENMP

\cp ${SOURCE_CROCO}/cppdefs.h cppdefs.h.MPI
sed 's/'undef\ \ \*RVTK_DEBUG'/'define\ RVTK_DEBUG'/' < cppdefs.h.MPI > cppdefs.h.MPI.tmp
mv cppdefs.h.MPI.tmp cppdefs.h.MPI
\cp ${SOURCE_CROCO}/param.h param.h.MPI


#List of test cases with only one points in one direction
LIST_2DV_X='GRAV_ADJ IGW INNERSHELF INTERNAL SHOREFACE SWASH THACKER TANK I_SOLITON KH_INST SANDBAR'
LIST_2DV_Y='OVERFLOW SHELFRONT'

source configure_file

Is2DV_X=0
Is2DV_Y=0
[ -n "$(echo $LIST_2DV_X |grep "${CONFIG_NAME}")" ] && Is2DV_X=1
[ -n "$(echo $LIST_2DV_Y |grep "${CONFIG_NAME}")" ] && Is2DV_Y=1

##############################################################################
#
#   FILL THE CPPDEFS.H 
# Title
echo TESTS OF $CONFIG_NAME
#
# 1- UNDEF ALL THE KEYS
#
echo 'undef '$LIST_KEY0
for par in SERIAL OPENMP MPI ; do 
    echo 'PARA=' $par
    for EXAMPLE in $LIST_KEY0
    do
	sed 's/'define\ \ \*$EXAMPLE'/'undef\ $EXAMPLE'/' < cppdefs.h.$par > cppdefs.h.$par.tmp
	\mv cppdefs.h.$par.tmp cppdefs.h.$par
    done
    
    # 2- DEFINE THE TYPE OF DEBUG TEST 
    sed 's/'undef\ \ \*$KEY_DEBUG'/'define\ $KEY_DEBUG'/' < cppdefs.h.$par > cppdefs.h.$par.tmp
    \mv cppdefs.h.$par.tmp cppdefs.h.$par
    
    
    # 3- DEFINE THE NAME OF THE CONFIG
    sed 's/'undef\ \*BENGUELA_LR'/'define\ $CONFIG_NAME'/' < cppdefs.h.$par > cppdefs.h.$par.tmp
    \mv cppdefs.h.$par.tmp cppdefs.h.$par
    
    # 4- DEFINE THE VARIOUS CPPKEYS
    #=4.1
    for EXAMPLE in $LIST_KEY_PHYS
    do
	sed 's/'undef\ \ \*$EXAMPLE'/'define\ $EXAMPLE'/' < cppdefs.h.$par > cppdefs.h.$par.tmp
	\mv cppdefs.h.$par.tmp cppdefs.h.$par
    done
    #==4.2
    for EXAMPLE in $LIST_KEY_NEST
    do
	echo $EXAMPLE
	sed 's/'undef\ \ \*$EXAMPLE'/'define\ $EXAMPLE'/' < cppdefs.h.$par > cppdefs.h.$par.tmp
	/bin/mv cppdefs.h.$par.tmp cppdefs.h.$par
	if [ "$LIST_KEY_NEST" = "AGRIF" ]; then
	    sed 's/'define\ \ \*AGRIF_2W'/'undef\ AGRIF_2W'/' < cppdefs.h.$par > cppdefs.h.$par.tmp
	    /bin/mv cppdefs.h.$par.tmp cppdefs.h.$par
	fi
    done
done

#=====================================================================================================
#=====================================================================================================
#####
echo 'FLAG OPENMP' ${FLAG_OPENMP}
echo 'FLAG MPI' ${FLAG_MPI}
####
#echo ' '
echo '==============================='
echo 'START TESTING ...             '
#echo '==============================='

SUCCESS=0
SUCCESS_COMP=0
SUCCESS_EXEC=0
if [ ! -f ${TEST_NAME}_steps ]; then 
    echo 'Y' > ${TEST_NAME}_steps
    echo 'Y' >> ${TEST_NAME}_steps
    echo 'Y' >> ${TEST_NAME}_steps
fi
#echo -e "   - Run Tests"> /dev/tty
#echo -e "   - Run Tests" > /dev/stdin

##############################################################################
# Serial runs
##############################################################################
echo "KEYS TESTED : "$LIST_KEY_PHYS
#echo ''
par1='SERIAL'
echo "SERIAL NPP=1 TEST $mytest"
[ -f check_file ] && /bin/rm check_file
sed 's/'NSUB_X=1,\ \ \*NSUB_E=NPP'/'NSUB_X=1,\ NSUB_E=1'/' < param.h.$par1 > param.h.$par1.tmp
sed 's/'NPP=1'/'NPP=1'/' < param.h.$par1 > param.h.$par1
\mv param.h.$par1.tmp param.h.$par1

[ -e  param.h.OK ] && \rm param.h.OK
[ -e  param.h.OK ] && \rm cppdefs.h.OK

\cp param.h.$par1 param.h.OK
\cp cppdefs.h.$par1 cppdefs.h.OK

Fqsub_serial
myreturn=$?

SUCCESS=$(($SUCCESS+$myreturn))
if [ "$myreturn" -eq 1 ]; then
    SUCCESS_COMP=$(($SUCCESS_COMP+1))
    sed -e '1c N' ${TEST_NAME}_steps > tmp.txt 
    \mv tmp.txt ${TEST_NAME}_steps 
    sed -e '2c ?' ${TEST_NAME}_steps > tmp.txt 
    \mv tmp.txt ${TEST_NAME}_steps 
fi  
if [ "$myreturn" -eq 2 ]; then
    SUCCESS_EXEC=$(($SUCCESS_EXEC+1))
    sed -e '2c N' ${TEST_NAME}_steps > tmp.txt 
    \mv tmp.txt ${TEST_NAME}_steps 
fi  

# 4- 
##############################################################################
# Openmp runs
##############################################################################
if [ ${FLAG_OPENMP} -eq 1 ]; then 
    
    par1='OPENMP'
    if [ $Is2DV_Y == 1 ]; then
	echo "OPEN-MP 1X2 NPP=2 TEST $mytest"
	sed 's/'NSUB_X=1,\ \ \*NSUB_E=NPP'/'NSUB_X=1,\ NSUB_E=2'/' < param.h.$par1 > param.h.$par1.tmp
	sed 's/'NPP=4'/'NPP=2'/' < param.h.$par1.tmp > param.h.$par1
    elif [ $Is2DV_X == 1 ]; then
	echo "OPEN-MP 2x1 NPP=2 TEST $mytest"
	sed 's/'NSUB_X=1,\ \ \*NSUB_E=NPP'/'NSUB_X=2,\ NSUB_E=1'/' < param.h.$par1 > param.h.$par1.tmp
	sed 's/'NPP=4'/'NPP=2'/' < param.h.$par1.tmp > param.h.$par1
    else
	echo "MPI ${NBPROCS_X}X${NBPROCS_Y} TEST $mytest"
	sed 's/'NP_XI=1,\ \ \*NP_ETA=4'/'NP_XI=${NBPROCS_X},\ NP_ETA=${NBPROCS_Y}'/' < param.h.$par1 > param.h.$par1.tmp
	sed 's/'NPP=4'/'NPP=4'/' < param.h.$par1.tmp > param.h.$par1
    fi
    \mv param.h.$par1.tmp param.h.$par1
    #
    sed 's/'undef\ \ \*${par1}'/'define\ ${par1}'/' < cppdefs.h.$par1 > cppdefs.h.$par1.tmp
    \mv cppdefs.h.$par1.tmp cppdefs.h.$par1
    #
    [ -e  param.h.OK ] && \rm param.h.OK
    [ -e  cppdefs.h.OK ] && \rm cppdefs.h.OK
    \cp param.h.$par1 param.h.OK
    \cp cppdefs.h.$par1 cppdefs.h.OK
    #
    Fqsub_openmp
    myreturn=$?

    SUCCESS=$(($SUCCESS+$myreturn))
    if [ "$myreturn" -eq 1 ]; then
	SUCCESS_COMP=$(($SUCCESS_COMP+1))
	sed -e '1c N' ${TEST_NAME}_steps > tmp.txt
	\mv tmp.txt ${TEST_NAME}_steps 
	sed -e '2c ?' ${TEST_NAME}_steps > tmp.txt 
	\mv tmp.txt ${TEST_NAME}_steps 
    fi  
    if [ "$myreturn" -eq 2 ]; then
	SUCCESS_EXEC=$(($SUCCESS_EXEC+1))
	sed -e '2c N' ${TEST_NAME}_steps > tmp.txt 
	\mv tmp.txt ${TEST_NAME}_steps 
    fi  

    #else
    
    #  sed -e '2c ?' ${TEST_NAME}_steps > tmp.txt 
    #  \mv tmp.txt ${TEST_NAME}_steps     
fi 

###############################################################################
# 4- RVTK_DEBUG_REG_DEV
##############################################################################
# Mpi runs
##############################################################################
if [ ${FLAG_MPI} -eq 1 ]; then 
    
    par1='MPI'
    if [ $Is2DV_Y == 1 ]; then
	echo "MPI 1X2 TEST $mytest"
	sed 's/'NP_XI=1,\ \ \*NP_ETA=4'/'NP_XI=1,\ NP_ETA=2'/' < param.h.$par1 > param.h.$par.tmp
	
    elif [ $Is2DV_X == 1 ]; then
	echo "MPI 2X1 TEST $mytest"
	sed 's/'NP_XI=1,\ \ \*NP_ETA=4'/'NP_XI=2,\ NP_ETA=1'/' < param.h.$par1 > param.h.$par.tmp
    else	
	echo "MPI ${NBPROCS_X}X${NBPROCS_Y} TEST $mytest"
	sed 's/'NP_XI=1,\ \ \*NP_ETA=4'/'NP_XI=${NBPROCS_X},\ NP_ETA=${NBPROCS_Y}'/' < param.h.$par1 > param.h.$par1.tmp
    fi
    \mv param.h.$par1.tmp param.h.$par1
    #
    sed 's/'undef\ \ \*$par1'/'define\ $par1'/' < cppdefs.h.$par1 > cppdefs.h.$par1.tmp
    \mv cppdefs.h.$par1.tmp cppdefs.h.$par1
    sed 's/'define\ \ \*MPI_NOLAND'/'undef\ MPI_NOLAND'/' < cppdefs.h.$par1 > cppdefs.h.$par1.tmp
    \mv cppdefs.h.$par1.tmp cppdefs.h.$par1
    #
    [ -e  param.h.OK ] && \rm param.h.OK
    [ -e  cppdefs.h.OK ] && \rm cppdefs.h.OK

    \cp param.h.$par1 param.h.OK
    \cp cppdefs.h.$par1 cppdefs.h.OK
    #
    Fqsub_mpi
    myreturn=$?
    
    SUCCESS=$(($SUCCESS+myreturn))
    
    if [ "$myreturn" -eq 1 ]; then
	SUCCESS_COMP=$(($SUCCESS_COMP+1))
	sed -e '1c N' ${TEST_NAME}_steps > tmp.txt 
	\mv tmp.txt ${TEST_NAME}_steps 
	sed -e '2c ?' ${TEST_NAME}_steps > tmp.txt 
	\mv tmp.txt ${TEST_NAME}_steps 
    fi  
    if [ "$myreturn" -eq 2 ]; then
	SUCCESS_EXEC=$(($SUCCESS_EXEC+1))
	sed -e '2c N' ${TEST_NAME}_steps > tmp.txt 
	\mv tmp.txt ${TEST_NAME}_steps 
    fi  
fi

if [  "$SUCCESS" -ne 0 ]; then
    sed -e '3c ?' ${TEST_NAME}_steps > tmp.txt ; \mv tmp.txt ${TEST_NAME}_steps
    echo  
    echo "SOMETHING WRONG HAPPENED"
    echo "EXITING ..."
    echo
    #  echo  > /dev/stdin
    #  echo -e "$(tput setaf 1 ; tput bold)SOMETHING WRONG HAPPENED WITH ${CONFIG_NAME} $(tput sgr0)" > /dev/stdin
    #  echo -e "$(tput setaf 1 ; tput bold)EXITING ...$(tput sgr0)"  > /dev/stdin
    #  echo  > /dev/stdin
    echo  | tee -a mylog.txt
    echo -e "${FMT_REDBLD}SOMETHING WRONG HAPPENED WITH ${CONFIG_NAME} ${FMT_ORD}" | tee -a mylog.txt
    echo -e "${FMT_REDBLD}EXITING ...${FMT_ORD}"  | tee -a mylog.txt 
    echo  | tee -a mylog.txt
    exit  1
fi

#########################################################################################################
# 5- Extract results
##############################################################################
#  runs
#echo ' '
echo "EXTRACTION $mytest"
Fextract_results $FLAG_MPI $FLAG_OPENMP
