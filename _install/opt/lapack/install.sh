#!/bin/bash
NAME=LAPACK
LIBNAME=liblapack.a

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

create_makeinc(){    
    local PLAT=$1
    cd $WRK_INSTALL
    DIR_TARGET=$PREFIX/$PLAT
    local BIN_TARGET=$DIR_TARGET/bin
    ETC_TARGET=$DIR_TARGET/etc
    local LIB_TARGET=$DIR_TARGET/lib
    local INC_TARGET=$DIR_TARGET/include
    OBJ_INSTALL=$SRC_INSTALL/obj_$PLAT
    echo "Creating directories:" >&2
    mkdir -pv $DIR_TARGET
    mkdir -pv $BIN_TARGET
    mkdir -pv $ETC_TARGET/modules/$LNAME
    mkdir -pv $LIB_TARGET
    mkdir -pv $INC_TARGET
    mkdir -pv $OBJ_INSTALL
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
    
    cat << EOF > make.inc
SHELL = /bin/sh
PLAT = $PLAT
OBJ_INSTALL=$OBJ_INSTALL
#  
#  Modify the FORTRAN and OPTS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader
#  and desired load options for your machine.
#
FORTRAN = $FC
OPTS     = $FFLAGS 
DRVOPTS  = $OPTS
NOOPT    = -O0
LOADER   = $FC
LOADOPTS = 

# Timer for the SECOND and DSECND routines
TIMER     = NONE
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= cvq
RANLIB   = ranlib
#
#  The location of BLAS library for linking the testing programs.
#  The target's machine-specific, optimized BLAS library should be
#  used whenever possible.
#
BLASLIB      = $HOME/BLAS/libblas.a
#
#  Location of the extended-precision BLAS (XBLAS) Fortran library
#  used for building and testing extended-precision routines.  The
#  relevant routines will be compiled and XBLAS will be linked only if
#  USEXBLAS is defined.
#
# USEXBLAS    = Yes
XBLASLIB     =
# XBLASLIB    = -lxblas
#
#  Names of generated libraries.
#
LAPACKLIB    = $LIB_TARGET/$LIBNAME
TMGLIB       = tmglib$PLAT.a
EIGSRCLIB    = eigsrc$PLAT.a
LINSRCLIB    = linsrc$PLAT.a

EOF

    echo "Copying init script for $UNAME" >&2
    cp -fv $BIN_INSTALL/configvars.sh $BIN_TARGET/configvars.sh
    cat <<EOF >> $BIN_TARGET/configvars.sh
add_library_to_system ${PREFIX}/${PLAT}
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

create_makeinc $PLAT
if [ ! -z $2 ];then
    if [ $2 == "clean" ];then
	make cleanall
	exit
    fi
fi
if [ -d $OBJ_INSTALL ];then
    rsync -av $OBJ_INSTALL/* $SRC_INSTALL/ 2>/dev/null
fi
echo "gUnzip SRC dir:";gunzip -r src
make
if [ $? == 0 ];then
    make clean
    mv -vf $WRK_INSTALL/make.inc $ETC_TARGET/make.inc.lapack
    echo "gZip SRC dir:" ; gzip -r src;gunzip src/Makefile.gz;gunzip src/VARIANTS/Makefile.gz
else
    echo "Error from Makefile. STOP here."
    exit 1
fi
rsync -av $WRK_INSTALL/share/man $PREFIX/$PLAT/

