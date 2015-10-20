#!/bin/bash
NAME=DMFT_TOOLS
LIBNAME=libdmftt.a

################### >>>> FUNCTIONS <<<<<< ###################
usage(){
    echo "usage:  $0 (-h,--help) [plat:intel,intel_debug,gnu,gnu_debug]"
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
VERSION=$(git describe --tags --always --long 2>/dev/null)
WRK_INSTALL=$WRKDIR/_install
if [ ! -d $WRK_INSTALL ];then echo "$0: can not find _install directory";exit;fi
BIN_INSTALL=$WRK_INSTALL/bin
ETC_INSTALL=$WRK_INSTALL/etc
ENVMOD_INSTALL=$ETC_INSTALL/environment_modules
SRC_INSTALL=$WRK_INSTALL/src
PREFIX=$(nparent_dir $WRK_INSTALL 4)

SFROOT=$PREFIX
if [ ! -d $SFROOT ];then echo "$0: can not find SciFor root directory";exit;fi
if [ ! -d $SFROOT/$PLAT ];then echo "$0: can not find $SFROOT/$PLAT directory";exit;fi


create_makeinc(){
    local PLAT=$1
    cd $WRK_INSTALL
    DIR_TARGET=$PREFIX/$PLAT
    local BIN_TARGET=$DIR_TARGET/bin
    ETC_TARGET=$DIR_TARGET/etc
    local LIB_TARGET=$DIR_TARGET/lib
    local INC_TARGET=$DIR_TARGET/include
    OBJ_INSTALL=$SRC_INSTALL/obj_$PLAT
    INC_INSTALL=$SRC_INSTALL/mod_$PLAT
    echo "Creating directories:" >&2
    mkdir -pv $DIR_TARGET
    mkdir -pv $BIN_TARGET
    mkdir -pv $ETC_TARGET/modules/$LNAME
    mkdir -pv $LIB_TARGET
    mkdir -pv $INC_TARGET
    mkdir -pv $OBJ_INSTALL
    mkdir -pv $INC_INSTALL
    case $PLAT in
	intel)
	    local FC=ifort
	    local FFLAGS='-O2 -ftz -static-intel'
	    local MOPT="-module "
	    ;;
	gnu)
	    local FC=gfortran
	    local FFLAGS='-O2 -funroll-all-loops -static'
	    local MOPT=-J
	    ;;
	intel_debug)
	    local FC=ifort
	    local FFLAGS='-p -O0 -g -debug -fpe0 -traceback -check all,noarg_temp_created -static-intel'
	    local MOPT="-module "
	    ;;
	gnu_debug)
	    FC=gfortran
	    FFLAGS='-O0 -p -g -Wall -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -static'
	    MOPT=-J
	    ;;
	ibm)
	    FC=xlf90
	    FFLAGS="-O1 -qarch=qp -qtune=qp"
	    MOPT="-qmoddir="
	    MOD_DIR=ibm_mods
	    OBJ_DIR=ibm_objs
	    ;;
	*)
	    usage
	    ;;
    esac
    
    cat << EOF > make.inc
FC=$FC
FFLAGS=$FFLAGS
MOPT=$MOPT
PLAT=$PLAT
LIB_DMFTT=$LIB_TARGET/libdmftt.a
INC_SCIFOR=$SFROOT/$PLAT/include
INC_TARGET=$INC_TARGET
OBJ_INSTALL=$OBJ_INSTALL
INC_INSTALL=$INC_INSTALL
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
make all
if [ $? == 0 ];then
    make clean
    mv -vf $WRK_INSTALL/make.inc $ETC_TARGET/make.inc.dmfttols
else
    echo "Error from Makefile. STOP here."
    exit 1
fi
