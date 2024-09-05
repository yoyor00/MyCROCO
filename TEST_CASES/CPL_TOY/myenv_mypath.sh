#!bin/bash
# This environnement works for Ubuntu/debian. Before starting
# you need to make sure that your system has the right compilers
# and libraries

# For CROCO, you will need:
#     - C compiler       (sudo apt-get install gcc )
#     - Fortran compiler (sudo apt-get install gfortran )
#     - Netcdf library   (sudo apt-get install libnetcdf-dev libnetcdff-dev )
#     - MPI library      (sudo apt-get install openmpi-bin libopenmpi-dev )
# For more informations please go here:
#    https://croco-ocean.gitlabpages.inria.fr/croco_doc/index.html 

# If you are coupling with WRF you can check this page:
#     https://www2.mmm.ucar.edu/wrf/OnLineTutorial/compilation_tutorial.php
# to see if you have the proper architecture
# the main libraries you need are ( in addtion to the above):
#     - libpng (sudo apt-get install libpng-dev)
#     - zlib   ( sudo apt-get instal zlib1g-dev)

# Make sure you also have nco tools!


#################################################################
################### FOR COMPILATIONS ############################
#################################################################
# -- for NETCDF
# -->
export NETCDF=$( nf-config --prefix )
export NETCDF_INC=$(nf-config --includedir)
export NETCDF_LIB=$( nf-config --flibs)

# -- Option for job launching
export MPI_LAUNCH=mpirun
export MPI_ext="--app"
export ncomod='nco'

#----------------------------------------------------------------
# Environment variables related to compilers
#----------------------------------------------------------------
export CC=gcc
export FC=gfortran
export F77=gfortran
export F90=gfortran
export MPIF90=mpif90
export MPICC=mpicc


#-------------
# WRF
#-------------
export NETCDF4=1
export NETCDF_classic=1
export WRFIO_NCD_LARGE_FILE_SUPPORT=1
#-------------
# WPS
#-------------
export JASPERLIB=""
export JASPERINC=""
#-------------
# WW3
#-------------
export NETCDF_CONFIG=${NETCDF}/bin/nf-config
export WWATCH3_NETCDF=NC4

#################################################################
####################### NEEDED PATHS ############################
#################################################################

#-------------------------------------------------------------------------------
# Machine settings
#-------------------------------------------------------------------------------
export MACHINE="Linux"

#-------------------------------------------------------------------------------
# Config paths
#-------------------------------------------------------------------------------
export CONFIG=BENGUELA_TOY
export CHOME="$CI_PROJECT_DIR/test/home/BENGUELA_TOY"
export CWORK="$CI_PROJECT_DIR/test/work/BENGUELA_TOY"

#-------------------------------------------------------------------------------
# Tools paths
#-------------------------------------------------------------------------------
export SCRIPTDIR=$CHOME/SCRIPTS_TOOLBOX

#-------------------------------------------------------------------------------
# Model source paths #Insert the full path ( do not use "~" for home )
#-------------------------------------------------------------------------------
export CPL="/oasis/compile_oa3-mct/"
export OCE="$CI_PROJECT_DIR/OCEAN"
export ATM="${HOME}/WRF"
export WAV="${HOME}/WW3/model"
export TOY="${CHOME}/TOY_IN"
export XIOS="${HOME}/XIOS"

#-------------------------------------------------------------------------------
# OASIS
#-------------------------------------------------------------------------------
export CPL_NAM_DIR="$CHOME/OASIS_IN"
export CPL_FILES_DIR="$CWORK/OASIS_FILES"

#-------------------------------------------------------------------------------
# OCE
#-------------------------------------------------------------------------------
export OCE_NAM_DIR="$CHOME/CROCO_IN"
export OCE_FILES_DIR="$CWORK/CROCO_FILES"
export OCE_FILES_ONLINEDIR=""

#-------------------------------------------------------------------------------
# ATM
#-------------------------------------------------------------------------------
export ATM_NAM_DIR="$CHOME/WRF_IN"
export ATM_FILES_DIR="$CWORK/WRF_FILES"

#-------------------------------------------------------------------------------
# WAV
#-------------------------------------------------------------------------------
export WAV_NAM_DIR="$CHOME/WW3_IN"
export WAV_FILES_DIR="$CWORK/WW3_FILES"

#-------------------------------------------------------------------------------
# TOY
#-------------------------------------------------------------------------------
export TOY_NAM_DIR="$CHOME/TOY_IN"

#------------------------------------------------------------------------------
# XIOS
#-------------------------------------------------------------------------------
export XIOS_NAM_DIR="$CHOME/XIOS_IN"

