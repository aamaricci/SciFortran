#!/bin/bash
NAME=PARPACK
LIBNAME=libparpack.a

################### >>>> FUNCTIONS <<<<<< ###################
usage(){
    echo "usage:  $0 [gnu,intel,gnu_debug,intel_debug]"
    exit
}


if [ -z $1 ] || [ $1 == "-h" ] || [ $1 == "--help"  ];then
    usage
fi

LIST_FC="gnu gnu_debug intel intel_debug"
[[ $LIST_FC =~ (^|[[:space:]])"$1"($|[[:space:]]) ]] || usage

nparent_dir(){
    local DIR=$1
    local N=$2
    for i in `seq 1 $N`;
    do 
	DIR=$(dirname $DIR)
    done
    echo $DIR
}

test_mpi(){
    local MPIROOT="$1"
    local TARGET="$2"
    local MPIF="$MPIROOT/include/mpif.h"
    if [ ! -e "$MPIF" ];then echo "Can not find the file $MPIF";exit;fi
    if [ ! -d "$TARGET/PARPACK/SRC/MPI" ];then echo "Can not find $TARGET/PARPACK/SRC/MPI";exit;fi
    if [ ! -d "$TARGET/PARPACK/UTIL/MPI" ];then echo "Can not find $TARGET/PARPACK/UTIL/MPI";exit;fi
    if [ ! -d "$TARGET/PARPACK/EXAMPLES/MPI" ];then echo "Can not find $TARGET/PARPACK/UTIL/MPI";exit;fi
    for DIR in SRC UTIL EXAMPLES
    do
	SFX=$(date  +%d_%m_%y_%Hh%M)
	FILE=$TARGET/PARPACK/$DIR/mpif.h
	if [ -e $FILE ];then
	    mv -vf $FILE ${FILE}.${SFX}
	fi
	cp -vf $MPIF $FILE
    done
}

################### >>>> PREAMBLE <<<<<< ###################
PLAT=$1
UNAME=`echo $NAME |tr [:lower:] [:upper:]`
LNAME=`echo $NAME |tr [:upper:] [:lower:]`
WRKDIR=$(pwd)
VERSION=$(git describe --tags 2>/dev/null)
WRK_INSTALL=$WRKDIR/_install
if [ ! -d $WRK_INSTALL ];then echo "$0: can not find _install directory";exit;fi
BIN_INSTALL=$WRK_INSTALL/bin
ETC_INSTALL=$WRK_INSTALL/etc
ENVMOD_INSTALL=$ETC_INSTALL/environment_modules
SRC_INSTALL=$WRK_INSTALL/SRC
PREFIX=$(nparent_dir $WRK_INSTALL 4)


which mpirun >/dev/null 2>&1
[[ $? == 0 ]] || usage
MPIROOT=$(nparent_dir $(which mpirun) 2)
test_mpi $MPIROOT $WRK_INSTALL

print_ARmake(){
    local PLAT=$1
    cd $WRK_INSTALL
    DIR_TARGET=$PREFIX/$PLAT
    ETC_TARGET=$DIR_TARGET/etc
    local BIN_TARGET=$DIR_TARGET/bin
    local LIB_TARGET=$DIR_TARGET/lib
    local INC_TARGET=$DIR_TARGET/include
    BLAS_INSTALL=$WRK_INSTALL/BLAS/blas_$PLAT
    LAPACK_INSTALL=$WRK_INSTALL/LAPACK/lapack_$PLAT
    UTIL_INSTALL=$WRK_INSTALL/UTIL/util_$PLAT
    OBJ_INSTALL=$SRC_INSTALL/obj_$PLAT
    PUTIL_INSTALL=$WRK_INSTALL/PARPACK/UTIL/MPI/util_$PLAT
    POBJ_INSTALL=$WRK_INSTALL/PARPACK/SRC/MPI/obj_$PLAT
    echo "Creating directories:" >&2
    mkdir -pv $DIR_TARGET
    mkdir -pv $BIN_TARGET
    mkdir -pv $ETC_TARGET/modules/$LNAME
    mkdir -pv $LIB_TARGET
    mkdir -pv $INC_TARGET
    mkdir -pv $OBJ_INSTALL $POBJ_INSTALL $BLAS_INSTALL $LAPACK_INSTALL $UTIL_INSTALL $PUTIL_INSTALL 
    case $PLAT in
	intel)
	    FC=ifort
	    FFLAGS='-O2 -ip -static-intel'
	    ;;
	gnu)
	    FC=gfortran
	    FFLAGS='-O2 -funroll-all-loops -static'
	    ;;
	intel_debug)
	    FC=ifort
	    FFLAGS='-p -O0 -g -debug -fpe0 -traceback -check all,noarg_temp_created -static-intel'
	    ;;
	gnu_debug)
	    FC=gfortran
	    FFLAGS='-O0 -p -g -Wall -fbacktrace -static'
	    ;;
	*)
	    usage
	    ;;
    esac
    
    cat << EOF > ARmake.inc
FC=$FC
FFLAGS=$FFLAGS
PLAT=$PLAT
OBJ_INSTALL=$OBJ_INSTALL
BLAS_INSTALL=$BLAS_INSTALL
LAPACK_INSTALL=$LAPACK_INSTALL
UTIL_INSTALL=$UTIL_INSTALL
PUTIL_INSTALL=$PUTIL_INSTALL
POBJ_INSTALL=$POBJ_INSTALL

home=$WRK_INSTALL
COMMLIB     = MPI

BLASdir      = \$(home)/BLAS
LAPACKdir    = \$(home)/LAPACK
UTILdir      = \$(home)/UTIL
SRCdir       = \$(home)/SRC
PSRCdir      = \$(home)/PARPACK/SRC/\$(COMMLIB)
PUTILdir     = \$(home)/PARPACK/UTIL/\$(COMMLIB)
DIRS   = \$(BLASdir) \$(LAPACKdir) \$(UTILdir) \$(SRCdir)


#  The name of the libraries to be created/linked to
ARPACKLIB  = $LIB_TARGET/libarpack.a
PARPACKLIB = $LIB_TARGET/libparpack.a
LAPACKLIB = 
BLASLIB = 
ALIBS =  \$(ARPACKLIB) \$(LAPACKLIB) \$(BLASLIB) 


# Libraries needed for Parallel ARPACK - MPI for SUN4
MPILIBS = #-L$MPIROOT/lib 
PLIBS = \$(PARPACKLIB) \$(ALIBS) \$(MPILIBS)


#  Make our own suffixes list.
#.SUFFIXES:
#.SUFFIXES:.f .o

#
#  Default command.
#
.DEFAULT:
	@\$(ECHO) "Unknown target \$@, try:  make help"

#
#  Command to build .o files from .f files.
#
.f.o:
	@\$(ECHO) Making \$@ from \$<
	@\$(FC) -c \$(FFLAGS) \$<

#
#  Various compilation programs and flags.
#  You need to make sure these are correct for your system.
LDFLAGS = 
CD	= cd
AR      = ar 
ARFLAGS  = cvq
CHMOD	= chmod
CHFLAGS	= -f
COMPRESS= compress

CP	= cp
CPP	 = /lib/cpp
CPPFLAGS =
ECHO	 = echo
LN	 = ln
LNFLAGS	 = -s
#MAKE	 = /bin/make
MKDIR	 = mkdir
MDFLAGS	 = -p
MV	 = mv
MVFLAGS	 = -f
RM	 = rm
RMFLAGS  = -f
SHELL	 = /bin/sh
TAR	 = tar
RANLIB   = ranlib


help:
	@\$(ECHO) "usage: make ?"

EOF
	
	echo "Copying init script for $UNAME" >&2
	cp -fv $BIN_INSTALL/configvars.sh $BIN_TARGET/configvars.sh
	cat <<EOF >> $BIN_TARGET/configvars.sh
add_library_to_system ${WRKDIR}/${PLAT}
EOF
	echo "" >&2
	echo "Generating environment module file for $UNAME" >&2
	cat <<EOF > $ETC_TARGET/modules/$LNAME/$PLAT
#%Modules
set	root	$WRKDIR
set	plat	$PLAT
set	version	"$VERSION ($PLAT)"
EOF
	cat $ENVMOD_INSTALL/module >> $ETC_TARGET/modules/$LNAME/$PLAT
	echo "" >&2
	echo "Compiling $UNAME library on platform $PLAT:">&2
	echo "" >&2
}


print_ARmake $PLAT
if [ ! -z $2 ];then
    if [ $2 == "clean" ];then
	make cleanall
	exit
    fi
fi
if [ -d $BLAS_INSTALL ];then
    rsync -av $BLAS_INSTALL/* $WRK_INSTALL/BLAS/
fi
if [ -d $LAPACK_INSTALL ];then
    rsync -av $LAPACK_INSTALL/* $WRK_INSTALL/LAPACK/
fi
if [ -d $UTIL_INSTALL ];then
    rsync -av $UTIL_INSTALL/* $WRK_INSTALL/UTIL/
fi
if [ -d $OBJ_INSTALL ];then
    rsync -av $OBJ_INSTALL/* $SRC_INSTALL/
fi
if [ -d $POBJ_INSTALL ];then
    rsync -av $POBJ_INSTALL/* $WRK_INSTALL/PARPACK/SRC/MPI/
fi
if [ -d $PUTIL_INSTALL ];then
    rsync -av $POBJ_INSTALL/* $WRK_INSTALL/PARPACK/UTIL/MPI/
fi
make all
if [ $? == 0 ];then
    make clean
    mv -vf $WRK_INSTALL/ARmake.inc $ETC_TARGET/make.inc.parpack
else
    echo "Error from Makefile. STOP here."
    exit 1
fi
