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

#
SCRDIR=$1
echo 'SRCDIR='$SCRDIR
#
# set source, compilation and run directories
#
#===
source CONFIGURE_GLOBAL
#===
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
FC=gfortran

#
# set MPI directories if needed
#
MPIF90="mpif90"
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
XIOS_ROOT_DIR=$HOME/xios-1.0
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
#####rm -rf $SCRDIR
mkdir $SCRDIR

#
# AGRIF sources directory
#
AGRIF_SRC=${ROOT_DIR}/AGRIF

#
# copy SOURCE code
#
/bin/cp -f ${SOURCE}/*.F90 $SCRDIR
/bin/cp -f ${SOURCE}/*.F   $SCRDIR
/bin/cp -f ${SOURCE}/*.h   $SCRDIR
/bin/cp -f ${SOURCE}/Make* $SCRDIR
/bin/cp -f ${SOURCE}/testkeys.F $SCRDIR
/bin/cp -f ${SOURCE}/jobcomp $SCRDIR
/bin/cp -f ${SOURCE}/amr.in $SCRDIR
/bin/cp -RLf ${AGRIF_SRC} $SCRDIR
/bin/cp -f ${ROOT_DIR}/XIOS/*.F $SCRDIR
/bin/cp -f ${ROOT_DIR}/PISCES/* $SCRDIR
/bin/cp -f ${ROOT_DIR}/PISCES/kRGB61* $RUNDIR
if [[ -e "namelist_pisces" ]] ; then
        echo "  file namelist_pisces exists in Run directory"
else
        /bin/cp -f ${SOURCE}/PISCES/namelist_pisces* $RUNDIR
        echo "  file namelist_pisces copied from source directory"
fi
#
# overwrite with local files
#
/bin/cp -f *.F90 $SCRDIR
/bin/cp -f *.F $SCRDIR
/bin/cp -f *.h $SCRDIR
/bin/cp -f Make* $SCRDIR
/bin/cp -f jobcomp $SCRDIR
#
#
# RVTK  files  DEBUG CPP KEYS
#
/bin/cp -f cppdefs_dev_cvtk.h ${SCRDIR}/cppdefs_dev.h
/bin/cp -f ${SCRDIR}/cppdefs.h.OK ${SCRDIR}/cppdefs.h
/bin/cp -f ${SCRDIR}/param.h.OK ${SCRDIR}/param.h
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
		FFLAGS1="-O3 -fno-alias -i4 -r8 -fp-model precise"
#                FFLAGS1="-O0 -g -i4 -r8 -traceback -check all -check bounds \
#                       -check uninit -CA -CB -CS -ftrapuv -fpe1"
		LDFLAGS1="$LDFLAGS1"
	elif [[ $FC == gfortran ]] ; then
		CPP1="cpp  -traditional -DLinux"
		CFT1=gfortran
		FFLAGS1="-O0 -fdefault-real-8 -fdefault-double-8 "
#		 FFLAGS1="-O0 -g -fdefault-real-8 -fdefault-double-8 -fbacktrace \
#			-fbounds-check -finit-real=nan -finit-integer=8888"
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
# determine if XIOS compilation is required
#
unset COMPILEXIOS
echo "Checking COMPILEXIOS..."
XIOS_ROOT_DIR=${CROCO_XIOS_ROOT_DIR-$XIOS_ROOT_DIR}
if $($CPP1 testkeys.F | grep -i -q xiosisdefined) ; then
        echo " => XIOS activated"
        COMPILEXIOS=TRUE
        LDFLAGS1="$LDFLAGS1 $XIOS_ROOT_DIR/lib/libxios.a  -lstdc++ -lnetcdff"
        CPPFLAGS1="$CPPFLAGS1 -I$XIOS_ROOT_DIR/inc"
        FFLAGS1="$FFLAGS1 -I$XIOS_ROOT_DIR/inc"
        ln -s $XIOS_ROOT_DIR/bin/xios_server.exe $RUNDIR/.
        $CPP1 -P -traditional -imacros cppdefs.h  ${ROOT_DIR}/XIOS/field_def.xml_full $RUNDIR/field_def.xml
        $CPP1 -P -traditional -imacros cppdefs.h  ${ROOT_DIR}/XIOS/domain_def.xml $RUNDIR/domain_def.xml
        $CPP1 -P -traditional -imacros cppdefs.h  ${ROOT_DIR}/XIOS/iodef.xml $RUNDIR/iodef.xml
fi

#
# determine if OASIS compilation is required
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
#
# compile the AGRIF librairy
#
	if [[ $COMPILEMPI ]] ; then
		$MAKE -C AGRIF FC="$CFT1" CPP="$CPP1" CPPFLAGS="-DAGRIF_MPI $MPIINC" FFLAGS="$FFLAGS1"
	else
		$MAKE -C AGRIF FC="$CFT1" CPP="$CPP1" FFLAGS="$FFLAGS1"
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
  
[ -f croco ] && mv croco $RUNDIR
#