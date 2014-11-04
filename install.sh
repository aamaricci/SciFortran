#!/bin/bash
NAME=DMFT_TOOLS
usage(){
    echo "usage:  $0 (-h,--help) [plat:intel,intel_debug,gnu,gnu_debug,ibm]"
    exit
}

if [ -z $1 ] || [ $1 == "-h" ] && [ $1=="--help"  ];then
    usage
fi

if [ $1 == "gnu" ] || [ $1 == "gnu_debug" ] || [ $1 == "ibm" ];then
    usage
fi

PLAT=$1
UNAME=`echo $NAME |tr [:lower:] [:upper:]`
LNAME=`echo $NAME |tr [:upper:] [:lower:]`
WRKDIR=$(pwd)
VERSION=$(git describe --tags)
WRK_INSTALL=$WRKDIR/_install
if [ ! -d $WRK_INSTALL ];then echo "$0: can not find _install directory";exit;fi

print_ARmake(){    
    cd $WRK_INSTALL
    local ROOT=$WRK_INSTALL
    local PLAT=$1
    local DIR_TARGET=$WRKDIR/$PLAT
    local BIN_TARGET=$DIR_TARGET/bin
    local ETC_TARGET=$DIR_TARGET/etc
    local LIB_TARGET=$DIR_TARGET/lib
    local INC_TARGET=$DIR_TARGET/include
    echo "Creating directories:" >&2
    mkdir -pv $DIR_TARGET
    mkdir -pv $BIN_TARGET
    mkdir -pv $ETC_TARGET;mkdir -pv $ETC_TARGET/modules
    mkdir -pv $LIB_TARGET
    mkdir -pv $INC_TARGET
    echo "" >&2
    case $PLAT in
	intel)
	    FC=ifort
	    FFLAGS="-O2 -static-intel"
	    MOPT="-module "
	    MOD_DIR=intel_mods
	    OBJ_DIR=intel_objs
	    ;;
	gnu)
	    FC=gfortran
	    FFLAGS="-O2 -static"
	    MOPT=-J
	    MOD_DIR=gnu_mods
	    OBJ_DIR=gnu_objs
	    ;;
	intel_debug)
	    FC=ifort
	    FFLAGS="-p -O0 -g -debug -fpe0 -traceback -check all,noarg_temp_created -static-intel"
	    MOPT="-module "
	    MOD_DIR=intel_debug_mods
	    OBJ_DIR=intel_debug_objs
	    ;;
	gnu_debug)
	    FC=gfortran
	    FFLAGS="-O0 -p -g -Wall -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -static"
	    MOPT=-J
	    MOD_DIR=gnu_debug_mods
	    OBJ_DIR=gnu_debug_objs
	    ;;
	ibm)
	    FC=xlf90
	    FFLAGS="-O2 -qarch=qp -qtune=qp"
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
OBJ_DIR=$ROOT/src/$OBJ_DIR
MOD_DIR=$ROOT/src/$MOD_DIR
LIB_DMFTT=$LIB_TARGET/libdmftt.a
MOD_DMFTT=$INC_TARGET
EOF


    echo "Copying init script for $UNAME" >&2
    cp -fv $WRK_INSTALL/bin/configvars.sh                         $BIN_TARGET/configvars.sh
    cat <<EOF >> $BIN_TARGET/configvars.sh
ROOT=$WRKDIR
PLAT=$PLAT
add_library_to_system $ROOT/$PLAT
EOF
    echo "" >&2
    echo "Generating environment module file for $UNAME" >&2
    cat <<EOF > $ETC_TARGET/modules/${LNAME}_$PLAT
#%Modules
set	root	$WRKDIR
set	plat	$PLAT
set	version	"$VERSION ($PLAT)"
EOF
    cat $WRK_INSTALL/etc/environment_modules/module >> $ETC_TARGET/modules/${LNAME}_$PLAT
    echo "" >&2
    echo "Compiling $UNAME library on platform $PLAT:">&2
    echo "" >&2

}


print_ARmake $PLAT
make all
if [ $? == 0 ];then
    make clean
    mv -vf make.inc $WRKDIR/$PLAT/
fi
cd $WRKDIR
