#!/bin/bash
NAME=SCIFOR
usage(){
    echo "usage:  $0 (-h,--help) [plat:intel,intel_debug,gnu,gnu_debug,ibm]"
    exit
}

if [ -z $1 ] || [ $1 == "-h" ] && [ $1=="--help"  ];then
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
    DIR_TARGET=$WRKDIR/$PLAT
    local BIN_TARGET=$DIR_TARGET/bin
    local ETC_TARGET=$DIR_TARGET/etc
    local LIB_TARGET=$DIR_TARGET/lib
    local INC_TARGET=$DIR_TARGET/include
    local LIB_SCIFOR=$LIB_TARGET/libscifor.a
    local MOD_SCIFOR=$INC_TARGET
    echo "Creating directories:" >&2
    mkdir -pv $DIR_TARGET
    mkdir -pv $BIN_TARGET
    mkdir -pv $ETC_TARGET
    mkdir -pv $ETC_TARGET/modules
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
LIB_SCIFOR=$LIB_SCIFOR
MOD_SCIFOR=$MOD_SCIFOR
EOF


    echo "Copying init script for $UNAME" >&2
    cp -fv $WRK_INSTALL/bin/scifor_completion.sh                  $BIN_TARGET/scifor_completion.sh
    cp -fv $WRK_INSTALL/bin/configvars.sh                         $BIN_TARGET/configvars.sh
    cat <<EOF >> $BIN_TARGET/configvars.sh
add_library_to_system ${WRKDIR}/${PLAT}
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


test_mkl(){
    local FILE="f.f90"
    echo "program test"     >  $FILE
    echo "end program test" >> $FILE
    local ARGS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm"
    echo "$FC $FILE $ARGS -o/dev/null"
    $FC $FILE $ARGS -o/dev/null
    local TEST=$?
    rm $FILE
    return $TEST
}



print_ARmake $PLAT
make all
if [ $? == 0 ];then
    make clean
    mv -vf make.inc $WRKDIR/$PLAT/
fi
cd $WRKDIR
CONFIGFILE=$WRKDIR/$PLAT/bin/configvars.sh
MODULEFILE=$WRKDIR/$PLAT/etc/modules/${LNAME}_$PLAT
mkdir -pv $HOME/.modules.d
mkdir -pv $HOME/.modules.d/applications
mkdir -pv $HOME/.modules.d/applications/$LNAME
cp -vf $MODULEFILE $HOME/.modules.d/applications/$LNAME/$PLAT


echo "" >&2
echo "TO USE LIBRARY $UNAME:" >&2
echo "A. add source $CONFIGFILE to your bash profile  (static)"   >&2
echo "     add the following line to your profile file (e.g. .bashrc)"      >&2
echo "     > echo \"source $CONFIGFILE\" "      >&2
echo "">&2
which modulecmd
if [ $? == 0 ];then
    echo "MODULES IS AVAILABLE IN YOUR SYSTEM SO YOU CAN USE THIS ALTERNATIVE METHOD:" >&2
    echo "B. ">&2
    echo "   load $MODULEFILE in your *environment module (dynamic).">&2
    echo "     add the following lines to your profile file (e.g. .bashrc)"     >&2
    echo "     > module use $HOME/.modules.d/applications">&2
    echo "     and then this to load the library: " >&2
    echo "     > module load $LNAME/$PLAT">&2
    echo "">&2
fi
