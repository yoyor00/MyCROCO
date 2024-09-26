export TERM="xterm-256color"
export CVTKHOME="/home/CVTK/test_repro"
export PERFRSTHOME="/home/CVTK/test_perfrst"
export DATADIR="/home/.datadir"
export PERFRSTDIR="/home/.perfrstdir"
export DATAREG="/data"
export CI_CROCO_PWD="/home"
export SOURCE_CROCO="/home/OCEAN"
export nest_position_reg="79 137 37 117 3 3 3 3"
export DATAVOR=""
export nest_position_vort=""
export DATAANA=""
export CI_FC="gfortran"
export CI_MPIF90="mpif90"
export CROCO_CI_MPIRUN="mpirun.openmpi --allow-run-as-root"
export CVTKWORK="/home/.datawork_gfortran"
export PERFRSTWORK="/home/.perfrstwork_gfortran"

export RUN_HOME=$CVTKHOME
export RUN_WORK=$CVTKWORK
export RUN_MK="mk_TESTALL"
export RUN_GATHER="gather_recap"

export RUN_NAME="REG"
export RUN_NAME_SCRIPT="reg"
export RUN_MK_NAME="reg"
export RUN_NAME_CONFIGURE="CONFIGURE_REG"

export CI_PROJECT_DIR="/home"

ulimit -s unlimited
mkdir -p $DATADIR
mkdir -p $CVTKWORK
mkdir -p $PERFRSTDIR
mkdir -p $PERFRSTWORK
cd $RUN_HOME/Scripts_${RUN_NAME_SCRIPT}
./create_link_master_${RUN_NAME_SCRIPT}.sh
cd -
cd $RUN_WORK/${RUN_NAME}
./${RUN_MK}.bash ${RUN_NAME_CONFIGURE} ${RUN_MK_NAME}
./${RUN_GATHER}.bash ${RUN_NAME} > /dev/null 2>&1
/bin/sh -c '! grep -i "Compilation failure" ${RUN_NAME}_gather_recap_*_git*' > /dev/null 2>&1
/bin/sh -c '! grep -i "Execution failure" ${RUN_NAME}_gather_recap_*_git*' > /dev/null 2>&1
/bin/sh -c '! grep -i "Parallel reproducibility failed" ${RUN_NAME}_gather_recap_*_git*' > /dev/null 2>&1  
 
