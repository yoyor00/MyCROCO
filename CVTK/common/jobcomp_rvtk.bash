#!/bin/bash
####################################################
#               COMPILATION JOB                    #
####################################################

# This script assumes default compilation options, to
# change those options : 
# it can either be  edited to add your own options
# or take into account the following 
# environment variables for compilation choices :
#
# CROCO_NETCDFLIB      : netcdf library
# CROCO_NETCDFINC      : netcdf include 
# CROCO_PRISM_ROOT_DIR : OASIS-MCT directory 
# CROCO_XIOS_ROOT_DIR  : XIOS directory
#
# CROCO_CFT1           : compiler
# CROCO_FFLAGS1        : compilation otpions
#
# Note that environment variables overwrite hard-coded
# options

#set -x
#
SCRDIR=$1
echo 'SRCDIR='$SCRDIR
#
# set source, compilation and run directories
#
SOURCE=$SOURCE_CROCO
echo 'SOURCE_CROCO='$SOURCE_CROCO
RUNDIR=`pwd`
ROOT_DIR=$SOURCE/..
#
# determine operating system
#
OS=`uname`
echo "OPERATING SYSTEM IS: $OS"

#
# compiler options
#
FC=$CI_FC

#
# set MPI directories if needed
#
MPIF90=$CI_MPIF90
MPILIB=""
MPIINC=""

#
# set NETCDF directories
#
#-----------------------------------------------------------
# Use : 
#-lnetcdf           : version netcdf-3.6.3                --
#-lnetcdff -lnetcdf : version netcdf-4.1.2                --
#-lnetcdff          : version netcdf-fortran-4.2-gfortran --
#-----------------------------------------------------------
#
#NETCDFLIB="-L/usr/local/lib -lnetcdf"
#NETCDFINC="-I/usr/local/include"
NETCDFLIB=$(nf-config --flibs)
NETCDFINC=-I$(nf-config --includedir)

#
# set OASIS-MCT (or OASIS3) directories if needed
#
PRISM_ROOT_DIR=../../../oasis3-mct/compile_oa3-mct

#
# set XIOS directory if needed
#
# if coupling with OASIS3-MCT is activated :
# => you need to use XIOS compiled with the "--use_oasis oasis3_mct" flag
#-----------------------------------------------------------
XIOS_ROOT_DIR=$HOME/xios

#
# END OF USER'S MODIFICATIONS
####################################################
#
# Use GNU Make command, else make
#
MAKE=gmake
which $MAKE > /dev/null 2>&1 || MAKE=make

#
# clean scratch area
#
rm -rf $SCRDIR
mkdir $SCRDIR

#
# AGRIF sources directory
#
AGRIF_SRC=${ROOT_DIR}/AGRIF

#
# copy SOURCE code
#

ls ${SOURCE}/*.F               > /dev/null  2>&1 && \cp ${SOURCE}/*.F   $SCRDIR
ls ${SOURCE}/*.F90             > /dev/null  2>&1 && \cp ${SOURCE}/*.F90 $SCRDIR
ls ${SOURCE}/*.h               > /dev/null  2>&1 && \cp ${SOURCE}/*.h   $SCRDIR
ls ${SOURCE}/Make*             > /dev/null  2>&1 && \cp ${SOURCE}/Make* $SCRDIR
ls ${SOURCE}/jobcomp           > /dev/null  2>&1 && \cp ${SOURCE}/jobcomp $SCRDIR
ls ${SOURCE}/amr.in            > /dev/null  2>&1 && \cp ${SOURCE}/amr.in $SCRDIR
ls ${AGRIF_SRC}                > /dev/null  2>&1 && \cp -r ${AGRIF_SRC} $SCRDIR
ls ${ROOT_DIR}/XIOS/*.F        > /dev/null  2>&1 && \cp ${ROOT_DIR}/XIOS/*.F $SCRDIR
ls ${ROOT_DIR}/PISCES/*        > /dev/null  2>&1 && \cp -r ${ROOT_DIR}/PISCES/* $SCRDIR
ls ${ROOT_DIR}/PISCES/SED/*    > /dev/null  2>&1 && \cp ${ROOT_DIR}/PISCES/SED/* $SCRDIR
ls ${ROOT_DIR}/PISCES/kRGB61*  > /dev/null  2>&1 && \cp ${ROOT_DIR}/PISCES/kRGB61* $RUNDIR
ls ${ROOT_DIR}/MUSTANG/*       > /dev/null  2>&1 && \cp -r ${ROOT_DIR}/MUSTANG/* $SCRDIR

if [[ -e "namelist_pisces_ref" ]] ; then
        echo "  file namelist_pisces exists in Run directory"
else
        \cp -f ${ROOT_DIR}/PISCES/namelist_pisces* $RUNDIR
        echo "  file namelist_pisces copied from source directory"
fi

if [[ -d MUSTANG_NAMELIST ]]; then
        echo "  Mustang namelist directory MUSTANG_NAMELIST exists"
else
        mkdir -p $RUNDIR/MUSTANG_NAMELIST
        \cp -rf ${ROOT_DIR}/MUSTANG/MUSTANG_NAMELIST/*txt $RUNDIR/MUSTANG_NAMELIST/.
        echo "  file para*txt copied from source directory"
fi

#
# overwrite with local files
#

ls *.F90   > /dev/null  2>&1 && \cp -f *.F90 $SCRDIR
ls *.F     > /dev/null  2>&1 && \cp -f *.F $SCRDIR
ls *.h     > /dev/null  2>&1 && \cp -f *.h $SCRDIR
ls *.h90   > /dev/null  2>&1 && \cp -f *.h90 $SCRDIR
ls Make*   > /dev/null  2>&1 && \cp -f Make* $SCRDIR
ls jobcomp > /dev/null  2>&1 && \cp -f jobcomp $SCRDIR
#
# RVTK  files  DEBUG CPP KEYS
#
/bin/cp -f cppdefs.h.OK ${SCRDIR}/cppdefs.h
/bin/cp -f param.h.OK ${SCRDIR}/param.h
#
# Change directory
#
cd $SCRDIR
#
# generates LDFLAGS1 according to users notifications
#
LDFLAGS1="${CROCO_NETCDFLIB-$NETCDFLIB}"
CPPFLAGS1="${CROCO_NETCDFINC-$NETCDFINC} -ICROCOFILES/AGRIF_INC"
#
# Set compilation options
#
if [[ $OS == Linux || $OS == Darwin ]] ; then           # ===== LINUX =====
	if [[ $FC == ifort || $FC == ifc ]] ; then
		CPP1="cpp -traditional -DLinux -DIfort"
		CFT1=ifort
                FFLAGS1="-O0 -mcmodel=medium -g -i4 -r8 -traceback -check all -check bounds \
                       -check uninit -CA -CB -CS -ftrapuv -fpe1"
		LDFLAGS1="$LDFLAGS1"
	elif [[ $FC == gfortran ]] ; then
		CPP1="cpp -traditional -DLinux"
		CFT1=gfortran
                FFLAGS1="-O0 -mcmodel=medium -g -fdefault-real-8 -fdefault-double-8 -std=legacy -fbacktrace \
			-fbounds-check -finit-real=zero -finit-integer=0"
		LDFLAGS1="$LDFLAGS1"
	fi
elif [[ $OS == CYGWIN_NT-10.0 ]] ; then  # ======== CYGWIN =======
        CPP1="cpp -traditional -DLinux"
        CFT1="gfortran"
        FFLAGS1="-O4 -fdefault-real-8 -fdefault-double-8 -march=native -mtune=native"
elif [[ $OS == AIX ]] ; then           # ===== IBM =====
	CPP1="cpp"
	CFT1="xlf95 -I$HOME/include/"
	MAKEAGRIF="Makefile.ibm.64"
	FFLAGS1="-q64 -qwarn64 -qfixed -qrealsize=8 -qintsize=8 -qhot \
			-qalias=noaryovrlp -qthreaded -O3 -qarch=pwr4 -qtune=pwr4 -qunroll=yes"
else
	echo "Unknown Operating System"
	exit
fi
#
# determine if AGRIF compilation is required
#
unset COMPILEAGRIF
echo "Checking COMPILEAGRIF..."
if $($CPP1 testkeys.F | grep -i -q agrifisdefined) ; then
	echo " => AGRIF activated"
	COMPILEAGRIF=TRUE
	FFLAGS1="$FFLAGS1 -IAGRIF"
	LDFLAGS1="-LAGRIF -lagrif $LDFLAGS1"
# we use the AGRIF Makedefs.generic definition
	cp -f Makedefs.generic.AGRIF Makedefs.generic
fi

#
# determine if MPI compilation is required
#
unset COMPILEMPI
echo "Checking COMPILEMPI..."
if $($CPP1 testkeys.F | grep -i -q mpiisdefined) ; then
	echo " => MPI activated"
	COMPILEMPI=TRUE
	LDFLAGS1="$LDFLAGS1 $MPILIB"
	CPPFLAGS1="$CPPFLAGS1 $MPIINC"
	FFLAGS1="$FFLAGS1 $MPIINC"
	CFT1="${MPIF90}"
fi

#
# Take environment variables for compiler and options
#
FFLAGS1=${CROCO_FFLAGS1-$FFLAGS1}
CFT1=${CROCO_CFT1-$CFT1}

#
# - Determine if XIOS librairies is required 
# - if it is the case :
#     => if XIOS compiled with oasis, add the OASIS inc. files and librairies
#     => pre-processing (using cpp) of the xml files required by XIOS 
#
unset COMPILEXIOS
echo "Checking COMPILEXIOS..."
XIOS_ROOT_DIR=${CROCO_XIOS_ROOT_DIR-$XIOS_ROOT_DIR}
if $($CPP1 testkeys.F | grep -i -q xiosisdefined) ; then
        echo " => XIOS activated"
        COMPILEXIOS=TRUE
        LDFLAGS1="$LDFLAGS1 $XIOS_ROOT_DIR/lib/libxios.a  -lstdc++ -lnetcdff -lnetcdf"
        CPPFLAGS1="$CPPFLAGS1 -I$XIOS_ROOT_DIR/inc"
        FFLAGS1="$FFLAGS1 -I$XIOS_ROOT_DIR/inc"
	
        ln -fs $XIOS_ROOT_DIR/bin/xios_server.exe $RUNDIR/.
fi

#
# determine if OASIS librairies are required
#
unset COMPILEOASIS
echo "Checking COMPILEOASIS..."
PRISM_ROOT_DIR=${CROCO_PRISM_ROOT_DIR-$PRISM_ROOT_DIR}
if $($CPP1 testkeys.F | grep -i -q oacplisdefined) ; then
    echo " => OASIS activated"
    CHAN=MPI1
    LIBPSMILE="${PRISM_ROOT_DIR}/lib/libpsmile.${CHAN}.a \
		${PRISM_ROOT_DIR}/lib/libmct.a  \
		${PRISM_ROOT_DIR}/lib/libmpeu.a \
		${PRISM_ROOT_DIR}/lib/libscrip.a"
    PSMILE_INCDIR="-I${PRISM_ROOT_DIR}/build/lib/psmile.${CHAN} \
		-I${PRISM_ROOT_DIR}/build/lib/mct"
    COMPILEOASIS=TRUE
    LDFLAGS1="$LDFLAGS1 $LIBPSMILE $NETCDFLIB"
    CPPFLAGS1="$CPPFLAGS1 ${PSMILE_INCDIR} $NETCDFINC"
    FFLAGS1="$FFLAGS1 ${PSMILE_INCDIR} $NETCDFINC"
fi
#
# prepare and compile the library
#
if [[ $COMPILEAGRIF ]] ; then
# Find the default C compiler
    CC1=$(echo -e 'dummy_target:\n\t@echo $(CC)' | $MAKE -f - dummy_target)
    CFLAGS1=$(echo -e 'dummy_target:\n\t@echo $(CFLAGS)' | $MAKE -f - dummy_target)
# Test if the C compiler is the GNU Compiler
# if True add '-fcommon' to CFLAGS
    if command -v "$CC1" >/dev/null && "$CC1" -v 2>&1 | grep -q "gcc version"; then
		echo "Using the GNU C compiler. Adding -fcommon to CFLAGS"
		CFLAGS1="$CFLAGS1 -fcommon"
	fi
#
# compile the AGRIF librairy
#
	if [[ $COMPILEMPI ]] ; then
		$MAKE -C AGRIF FC="$CFT1" CPP="$CPP1" CPPFLAGS="-DAGRIF_MPI $MPIINC" FFLAGS="$FFLAGS1" CFLAGS="$CFLAGS1"
	else
		$MAKE -C AGRIF FC="$CFT1" CPP="$CPP1" FFLAGS="$FFLAGS1" CFLAGS="$CFLAGS1"
	fi
	if [[ $OS == Darwin ]] ; then          # DARWIN
# run RANLIB on Darwin system
		ranlib AGRIF/libagrif.a
	fi
#
	mkdir CROCOFILES
	mkdir -p CROCOFILES/AGRIF_MODELFILES
	mkdir -p CROCOFILES/AGRIF_INC
	$CPP1 amr.in | grep -v -e ! -e '#' -e % -e '*' > CROCOFILES/amr.scrum
	mv AGRIF/conv CROCOFILES/.
	for i in *.h *.h90 ; do
		echo $i
		cat cppdefs.h $i | cpp -P | grep -v -e ! -e '#' -e % -e '*' > CROCOFILES/$i
	done
	mv -f CROCOFILES/private_scratch_AMR.h CROCOFILES/private_scratch.h
fi

#
# determine if OPENMP compilation is needed
#
unset COMPILEOMP
echo "Checking COMPILEOMP..."
if $($CPP1 testkeys.F | grep -i -q openmp) ; then
	COMPILEOMP=TRUE
	if [[ $OS == Linux || $OS == Darwin ]] ; then 
		if [[ $FC == gfortran ]] ; then
			FFLAGS1="$FFLAGS1 -fopenmp"
		elif [[ $FC == ifort || $FC == ifc ]] ; then
			FFLAGS1="$FFLAGS1 -openmp"
		else
			FFLAGS1="$FFLAGS1 -openmp"
		fi
	elif [[ $OS == CYGWIN_NT-10.0 ]] ; then
        FFLAGS1=="$FFLAGS1 -fopenmp"
	elif [[ $OS == AIX ]] ; then
		FFLAGS1="$FFLAGS1 -qsmp=omp"
		CFT1="xlf95_r"
	fi
fi

#
# rewrite Makedefs according to previous flags
# with openmp flags if needed
#
rm -f Makedefs
echo 's?$(FFLAGS1)?'$FFLAGS1'?g' > flags.tmp
echo 's?$(LDFLAGS1)?'$LDFLAGS1'?g' >> flags.tmp
echo 's?$(CPP1)?'$CPP1'?g' >> flags.tmp
echo 's?$(CFT1)?'$CFT1'?g' >> flags.tmp
echo 's?$(CPPFLAGS1)?'$CPPFLAGS1'?g' >> flags.tmp
sed -f flags.tmp Makedefs.generic > Makedefs
rm -f flags.tmp

#
# compile croco
#
$MAKE depend
$MAKE
  
[[ -f croco  ]] && mv croco $RUNDIR
#[[ -f partit ]] && mv partit $RUNDIR
#[[ -f ncjoin ]] && mv ncjoin  $RUNDIR
#
