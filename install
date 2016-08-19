#!/bin/bash
NAME=SCIFOR
UNAME=`echo $NAME |tr [:lower:] [:upper:]`
LNAME=`echo $NAME |tr [:upper:] [:lower:]`
LIBNAME=lib$LNAME.a
LOG=install.log
>$LOG
exec >  >(tee -a $LOG)
exec 2> >(tee -a $LOG >&2)

echo ""
echo "Invoked as:"
echo "$0 $*"
echo ""

USER=$(id -un)
GUSER=$(id -gn)


#>>> GET THE ACTUAL DIRECTORY
HERE=$(pwd)

#>>> USAGE FUNCTION
usage(){
    echo ""
    echo "usage:"
    echo ""
    echo "$0  --plat=FC_PLAT --prefix=PREFIX_DIR  [ -d,--debug   -h,--help   -m,--no-mpi   -q,--quiet ]"
    echo ""
    echo "mandatory arguments:" 
    echo "    --plat      : specifies the actual platform/compiler to use [default: intel] (valid: intel, gnu)"
    echo "    --prefix    : specifies the target directory [default: ~/opt/scifor]"
    echo ""    
    echo "optional arguments:"
    echo "    -d,--debug  : debug flag"
    echo "    -h,--help   : this help"
    echo "    -m,--no-mpi : do not use MPI compiler: produce a serial library"
    echo "    -q,--quiet  : assume Y to all questions."
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
params="$(getopt -n "$0" --options dhmq --longoptions plat:,prefix:,debug,help,no-mpi,quiet -- "$@")"
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


#>>> SET SOME DEFAULTS VARIABLES AND OTHER ONES.  0=True; 1=False
DEBUG=1				# Assume DEBUG=False
QUIET=1				# Assume QUIET=False
IDEFAULT=0			# Assume default installation (True)
WMPI=0				# Assume with MPI (WMPI=True)
OPT=""
#DEFAULT PREFIX
PREFIX=$HOME/opt/scifor
#DEFAULT PLATFORM
PLAT=gnu
VERSION=$(git describe --tags 2>/dev/null)
WRK_INSTALL=$(pwd)
BIN_INSTALL=$WRK_INSTALL/bin
ETC_INSTALL=$WRK_INSTALL/etc
OPT_INSTALL=$WRK_INSTALL/opt
ENVMOD_INSTALL=$ETC_INSTALL/environment_modules
SRC_INSTALL=$WRK_INSTALL/src


#>>> A COUPLE OF LISTS OF ALLOWED OPTIONS
LIST_FC="gnu intel"


#>>> GO THROUGH THE INPUT ARGUMENTS. FOR EACH ONE IF REQUIRED TAKE ACTION BY SETTING VARIABLES.
while true
do
    case $1 in
	--plat)
	    IDEFAULT=1
	    PLAT=$2
	    shift 2
	    [[ ! $LIST_FC =~ (^|[[:space:]])"$PLAT"($|[[:space:]]) ]] && {
		echo "Incorrect Fortran PLAT: $PLAT";
		echo " available values are: $LIST_FC"
		exit 1
	    }
	    ;;
	--prefix)
	    IDEFAULT=1
	    PREFIX=$2;
	    shift 2
	    ;;
	-d|--debug) DEBUG=0;shift ;;
        -h|--help) usage ;;
	-m|--no-mpi) WMPI=1;shift ;;
	-q|--quiet) QUIET=0;shift ;;
        --) shift; break ;;
        *) usage ;;
    esac
done

#>>> CHECK THAT THE MANDATORY OPTIONS ARE PRESENT:
[[ $IDEFAULT == 0 ]] && {
    LIST_ARGS+=" --prefix=$PREFIX --plat=$PLAT"
    echo "Using default configuaration --prefix=$PREFIX, --plat=$PLAT"; sleep 0.5;echo ""
}

[[ ! -z $PLAT ]] && [[ ! -z $PREFIX ]] || usage


case $PLAT in
    gnu) _FC_=gfortran ;;
    intel) _FC_=ifort ;;
esac


#RENAME WITH DEBUG IF NECESSARY 
[[ $DEBUG == 0 ]] && PLAT=${PLAT}_debug
[[ $WMPI == 1 ]] && PLAT=${PLAT}_nompi


#>>> SET STANDARD NAMES FOR THE TARGET DIRECTORY
DIR_TARGET=$PREFIX/$PLAT
BIN_TARGET=$DIR_TARGET/bin
ETC_TARGET=$DIR_TARGET/etc
LIB_TARGET=$DIR_TARGET/lib
INC_TARGET=$DIR_TARGET/include
DIR_TARGET_W=0
if [ -d $PREFIX ];then
    TEST_W_DIR=$PREFIX
else
    TEST_W_DIR=$(nparent_dir $PREFIX 1)
fi
if [ ! -w $TEST_W_DIR ];then
    DIR_TARGET_W=1
    echo "Can not create $DIR_TARGET: $TEST_W_DIR has no write access"
    sleep 0.5
    echo "Try to grant root privileges to create $DIR_TARGET for $USER:$GROUP"
    sudo -v
fi


#>>> SETUP MPI
MPIFC=mpif90
if [ $WMPI == 0 ];then
    which $MPIFC >/dev/null 2>&1
    CHECK_MPIFC=$?
    if [ $CHECK_MPIFC == 0 ];then
	echo "Using the following MPI compiler:"
	mpif90 -show
	sleep 0.5
	echo "On platform: $PLAT"
	sleep 0.5
	_MPIFC_=$(mpif90 -show | awk '{print $1}') # get the MPI f90 compiler
	if [ "$_MPIFC_" != "$_FC_" ];then
	    echo "Possible mismatch between MPI:+$_MPIFC_ and Platform $PLAT:+$_FC_"
	    echo "Check your MPI installation... exiting"
	    sleep 2
	    exit 1
	fi
    else
	echo "Can not find the $MPIFC compiler in the system."
	sleep 0.5
	echo "Trying installing with the --no-mpi option"
	sleep 2
	exit 1 
    fi
    CPP_FLAG='MPI'
else
    CPP_FLAG=''
fi


#>>> CREATE THE MAKE.INC IN THE FOLLOWING FUNCTION
create_makeinc(){
    local PLAT=$1
    cd $WRK_INSTALL
    case $PLAT in
	intel|intel_nompi)
	    local FC=ifort
	    local FCOPT="-O2 -ftz -static-intel -fpp -D_$CPP_FLAG "
	    local MOPT="-module "
	    ;;
	gnu|gnu_nompi)
	    local FC=gfortran
	    local FCOPT="-O2 -funroll-all-loops -static -cpp -D_$CPP_FLAG "
	    local MOPT=-J
	    ;;
	intel_debug|intel_debug_nompi)
	    local FC=ifort
	    local FCOPT="-p -O0 -g -debug -fpe0 -traceback -check all,noarg_temp_created -static-intel -fpp -D_$CPP_FLAG "
	    local MOPT="-module "
	    ;;
	gnu_debug|gnu_debug_nompi)
	    local FC=gfortran
	    local FCOPT="-O0 -p -g -Wall -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -static -cpp -D_$CPP_FLAG "
	    local MOPT=-J
	    ;;
	# ibm)
	#     FC=xlf90
	#     FCOPT="-O2 -qarch=qp -qtune=qp"
	#     MOPT="-qmoddir="
	#     ;;
	*)
	    usage
	    ;;
    esac
    
    if [ $WMPI == 0 ];then
	FC=$MPIFC		# set Fortran Compiler to the MPI Compiler (we checked it exists)
    fi

    cat << EOF > make.inc
FC=$FC
MPIFC=$MPIFC
FCOPT=$FCOPT
MOPT=$MOPT
PLAT=$PLAT
INC_TARGET=$INC_TARGET
LIB_SCIFOR=$LIB_TARGET/$LIBNAME
EOF
}

#>>> GET THE ACTUAL DIRECTORY
HERE=$(pwd)

create_makeinc $PLAT
sleep 2


if [ $QUIET == 1 ];then		# if NO quiet
    _DIR=Y
    echo -n "Installing in $DIR_TARGET. Continue [Y/n]: "
    read _DIR;
    _DIR=`echo $_DIR |tr [:lower:] [:upper:]`
    [[ $_DIR == Y ]] || exit 1
else
    echo "Installing SCIFOR in $DIR_TARGET (quiet mode): "
    sleep 2
fi




# >>> CREATE THE DIRECTORY HIERARCHY:
echo "Creating directories:"
sleep 1
if [ $DIR_TARGET_W -eq 0 ];then
    echo "mkdir -pv $DIR_TARGET"
    mkdir -pv $DIR_TARGET
else
    echo "sudo mkdir -pv $DIR_TARGET && sudo chown $USER:$GROUP $DIR_TARGET"
    sudo mkdir -pv $DIR_TARGET && sudo chown $USER:$GROUP $DIR_TARGET
fi
mkdir -pv $BIN_TARGET
mkdir -pv $ETC_TARGET/modules/$LNAME
mkdir -pv $LIB_TARGET
mkdir -pv $INC_TARGET
sleep 2


cd $HERE


echo "Copying init script for $UNAME" 
cp -fv $BIN_INSTALL/scifor_completion.sh  $BIN_TARGET/scifor_completion.sh
cp -fv $BIN_INSTALL/configvars.sh $BIN_TARGET/configvars.sh
cat <<EOF >> $BIN_TARGET/configvars.sh
add_library_to_system ${PREFIX}/${PLAT}
EOF
echo "" 
sleep 1


echo "Generating environment module file for $UNAME" 
cat <<EOF > $ETC_TARGET/modules/$LNAME/$PLAT
#%Modules
set	root	$PREFIX
set	plat	$PLAT
set	version	"$VERSION ($PLAT)"
set     compiler $FC
EOF
cat $ENVMOD_INSTALL/module >> $ETC_TARGET/modules/$LNAME/$PLAT
echo "" 
sleep 1

echo "Compiling $UNAME library on platform $PLAT:"
echo "" 
sleep 1

#CREATING THE SCIFOR LIBRARY:
#rm -fv $LIB_TARGET/$LIBNAME
if [ $WMPI == 0 ];then
    make parallel
else
    make serial
fi

if [ $? == 0 ];then
    make clean
    mv -vf $WRK_INSTALL/make.inc $ETC_TARGET/make.inc.scifor
    mv -vf $WRK_INSTALL/$LOG $DIR_TARGET/
else
    echo "Error from Makefile. STOP here."
    exit 1
fi
sleep 1


#LAST TOUCH COPY THE CONFIGVARS AND CREATE THE USER MODULES FILE. PRINT USAGE DETAILS.
CONFIGFILE=$PREFIX/$PLAT/bin/configvars.sh
MODULEFILE=$PREFIX/$PLAT/etc/modules/$LNAME/$PLAT
mkdir -pv $HOME/.modules.d/other/$LNAME
cp -vf $MODULEFILE $HOME/.modules.d/other/$LNAME/$PLAT

echo ""
echo "" 
echo "USAGE:" 
echo "" 
echo " To add SciFor (and other libraries) to your system:" 
echo "   $ source $CONFIGFILE" 
echo " (or copy this line into your bash profile [e.g. .bashrc])" 
echo ""
module avail >/dev/null 2>&1 
if [ $? == 0 ];then
echo " or load the $UNAME modules:" 
echo "   $ module use $HOME/.modules.d/other" 
echo "   $ module load $LNAME/$PLAT" 
echo "(or copy these lines into your bash profile [e.g. .bashrc])" 
echo ""
fi
echo ""
echo "Enjoy... (for info: adriano.amaricciATgmail.com)"
echo ""
cd $HERE
