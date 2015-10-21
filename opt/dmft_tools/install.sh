#!/bin/bash
NAME=DMFTT
UNAME=`echo $NAME |tr [:lower:] [:upper:]`
LNAME=`echo $NAME |tr [:upper:] [:lower:]`
LIBNAME=lib$LNAME.a
LOG=install.log
>$LOG
exec >  >(tee -a $LOG)
exec 2> >(tee -a $LOG >&2)


#>>> USAGE FUNCTION
usage(){
    echo ""
    echo "usage:"
    echo ""
    echo "$0  -p,--plat=FC_PLAT  [ --prefix=PREFIX_DIR  -c,--clean -d,--debug  -h,--help ]"
    echo ""
    echo "    -p,--plat   : specifies the actual platform/compiler to use [intel,gnu]"
    echo "    --prefix    : specifies the target directory [default: FC_PLAT]"
    echo "    -q,--quiet  : assume Y to all questions."
    echo "    -c,--clean  : clean out the former compilation."
    echo "    -d,--debug  : debug flag"
    echo "    -h,--help   : this help"
    echo ""
    exit
}



#>>> GET Nth TIMES PARENT DIRECTORY
nparent_dir(){
    local DIR=$1
    local N=$2
    for i in `seq 1 $N`;
    do 
	DIR=$(dirname $DIR)
    done
    echo $DIR
}

#>>> GET THE ENTIRE LIST OF ARGUMENTS PASSED TO STDIN
LIST_ARGS=$*

#>>> GET LONG & SHORT OPTIONS
params="$(getopt -n "$0" --options p:o:qcwdh --longoptions plat:,prefix:,opt-lib:,quiet,clean,wdmftt,debug,help -- "$@")"
if [ $? -ne 0 ];then
    usage
fi
eval set -- "$params"
unset params

#>>> CHECK THE NUMBER OF ARGUMENTS. IF NONE ARE PASSED, PRINT HELP AND EXIT.
NUMARGS=$#
if [ $NUMARGS -eq 0 ]; then
    usage
fi



#>>> SET SOME DEFAULTS VARIABLES AND OTHER ONES
WPLAT=1
DEBUG=1
CLEAN=1
WRK_INSTALL=$(pwd)
PREFIX=$(nparent_dir $WRK_INSTALL 2)
BIN_INSTALL=$WRK_INSTALL/bin
ETC_INSTALL=$WRK_INSTALL/etc
OPT_INSTALL=$WRK_INSTALL/opt
ENVMOD_INSTALL=$ETC_INSTALL/environment_modules
SRC_INSTALL=$WRK_INSTALL/src

#>>> THE LISTS OF ALLOWED PLAT
LIST_FC="gnu intel"


#>>> GO THROUGH THE INPUT ARGUMENTS. FOR EACH ONE IF REQUIRED TAKE ACTION BY SETTING VARIABLES.
while true
do
    case $1 in
	-p|--plat)
	    WPLAT=0
	    PLAT=$2
	    shift 2
	    [[ ! $LIST_FC =~ (^|[[:space:]])"$PLAT"($|[[:space:]]) ]] && {
		echo "Incorrect Fortran PLAT: $PLAT";
		echo " available values are: $LIST_FC"
		exit 1
	    }
	    ;;
	--prefix)
	    PREFIX=$2;
	    shift 2
	    ;;
	-c|--clean) CLEAN=0;shift ;;
	-d|--debug) DEBUG=0;shift ;;
        -h|--help) usage ;;
	-o|--opt-lib) shift 2;;
	-q|--quiet) shift ;;
	-w|--wdmftt) shift;;
        --) shift; break ;;
        *) usage ;;
    esac
done

#>>> CHECK THAT THE MANDATORY OPTION -p,-plat IS PRESENT:
[[ $WPLAT == 0 ]] || usage


#RENAME WITH DEBUG IF NECESSARY 
if [ $DEBUG == 0 ];then 
    PLAT=${PLAT}_debug;
fi

#>>> SET STANDARD NAMES FOR THE TARGET DIRECTORY
DIR_TARGET=$PREFIX/$PLAT
BIN_TARGET=$DIR_TARGET/bin
ETC_TARGET=$DIR_TARGET/etc
LIB_TARGET=$DIR_TARGET/lib
INC_TARGET=$DIR_TARGET/include
echo "Installing in $DIR_TARGET."
sleep 2


SFROOT=$PREFIX
if [ ! -d $SFROOT ];then echo "$0: can not find SciFor root directory";exit;fi
if [ ! -d $SFROOT/$PLAT ];then echo "$0: can not find $SFROOT/$PLAT directory";exit;fi


create_makeinc(){
    local PLAT=$1
    cd $WRK_INSTALL
    OBJ_INSTALL=$SRC_INSTALL/obj_$PLAT
    INC_INSTALL=$SRC_INSTALL/mod_$PLAT
    echo "Creating directories:" 
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

    echo "Copying init script for $UNAME" 
    cp -fv $BIN_INSTALL/configvars.sh $BIN_TARGET/configvars.sh
    cat <<EOF >> $BIN_TARGET/configvars.sh
add_library_to_system ${PREFIX}/${PLAT}
EOF
    echo "" 
    echo "Generating environment module file for $UNAME" 
    cat <<EOF > $ETC_TARGET/modules/$LNAME/$PLAT
#%Modules
set	root	$PREFIX
set	plat	$PLAT
set	version	"($PLAT)"
EOF
    cat $ENVMOD_INSTALL/module >> $ETC_TARGET/modules/$LNAME/$PLAT
    echo "" 
    echo "Compiling $UNAME library on platform $PLAT:"
    echo "" 
}

create_makeinc $PLAT
if [ $CLEAN == 0 ];then
    make cleanall
    exit 0
fi
make all
if [ $? == 0 ];then
    make clean
    mv -vf $WRK_INSTALL/make.inc $ETC_TARGET/make.inc.dmfttols
else
    echo "Error from Makefile. STOP here."
    exit 1
fi
