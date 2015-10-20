#!/bin/bash
NAME=ARPACK
LIBNAME=libarpack.a

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
SRC_INSTALL=$WRK_INSTALL/src
PREFIX=$(nparent_dir $WRK_INSTALL 4)

print_ARmake(){    
    local PLAT=$1
    cd $WRK_INSTALL
    DIR_TARGET=$PREFIX/$PLAT
    local BIN_TARGET=$DIR_TARGET/bin
    local LIB_TARGET=$DIR_TARGET/lib
    local INC_TARGET=$DIR_TARGET/include
    ETC_TARGET=$DIR_TARGET/etc
    BLAS_INSTALL=$WRK_INSTALL/blas/blas_$PLAT
    LAPACK_INSTALL=$WRK_INSTALL/lapack/lapack_$PLAT
    UTIL_INSTALL=$WRK_INSTALL/util/util_$PLAT
    OBJ_INSTALL=$SRC_INSTALL/obj_$PLAT
    echo "Creating directories:" >&2
    mkdir -pv $DIR_TARGET
    mkdir -pv $BIN_TARGET
    mkdir -pv $ETC_TARGET/modules/$LNAME
    mkdir -pv $LIB_TARGET
    mkdir -pv $INC_TARGET
    mkdir -pv $OBJ_INSTALL $BLAS_INSTALL $LAPACK_INSTALL $UTIL_INSTALL
    case $PLAT in
	intel)
	    local FC=ifort
	    local FFLAGS='-O2 -ftz -static-intel'
	    local MOPT=-module 
	    ;;
	gnu)
	    local FC=gfortran
	    local FFLAGS='-O2 -funroll-all-loops -static'
	    local MOPT=-J
	    ;;
	intel_debug)
	    local FC=ifort
	    local FFLAGS='-p -O0 -g -debug -fpe0 -traceback -check all,noarg_temp_created -static-intel'
	    local MOPT=-module 
	    ;;
	gnu_debug)
	    local FC=gfortran
	    local FFLAGS='-O0 -p -g -Wall -fbacktrace -static'
	    local MOPT=-J
	    ;;
	*)
	    usage
	    ;;
    esac
    
    rm -vf ARmake.inc
    cat << EOF >> ARmake.inc
FC=$FC
FFLAGS=$FFLAGS
PLAT=$PLAT
OBJ_INSTALL=$OBJ_INSTALL
BLAS_INSTALL=$BLAS_INSTALL
LAPACK_INSTALL=$LAPACK_INSTALL
UTIL_INSTALL=$UTIL_INSTALL

home=$WRK_INSTALL

BLASdir      = \$(home)/blas
LAPACKdir    = \$(home)/lapack
UTILdir      = \$(home)/util
SRCdir       = \$(home)/src
DIRS   = \$(BLASdir) \$(LAPACKdir) \$(UTILdir) \$(SRCdir)


#  The name of the libraries to be created/linked to
ARPACKLIB  = $LIB_TARGET/$LIBNAME
LAPACKLIB = 
BLASLIB = 
ALIBS =  \$(ARPACKLIB) \$(LAPACKLIB) \$(BLASLIB) 


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
set	root	$PREFIX
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
    rsync -av $BLAS_INSTALL/* $WRK_INSTALL/blas/
fi
if [ -d $LAPACK_INSTALL ];then
    rsync -av $LAPACK_INSTALL/* $WRK_INSTALL/lapack/
fi
if [ -d $UTIL_INSTALL ];then
    rsync -av $UTIL_INSTALL/* $WRK_INSTALL/util/
fi
if [ -d $OBJ_INSTALL ];then
    rsync -av $OBJ_INSTALL/* $SRC_INSTALL/
fi
make all
if [ $? == 0 ];then
    make clean
    mv -vf $WRK_INSTALL/ARmake.inc $ETC_TARGET/make.inc.arpack
else
    echo "Error from Makefile. STOP here."
    exit 1
fi
