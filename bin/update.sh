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
UPDATE=1			# Assume UPDATE=False (full installation default)
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
	-u|--update) UPDATE=0;shift ;;
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
if [ ! -d $DIR_TARGET ];then
    echo "Directory: $DIR_TARGET does not exist. UPDATE is impossible. Check or install from scratch."
    exit 1
else
    TEST_W_DIR=$DIR_TARGET
fi
if [ ! -w $TEST_W_DIR ];then
    DIR_TARGET_W=1
    echo "$TEST_W_DIR has no write access"
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



#>>> COPY BACK THE MAKE.INC
cp -vf $ETC_TARGET/make.inc.scifor $WRK_INSTALL/make.inc 
sleep 2


if [ $QUIET == 1 ];then		# if NO quiet
    _DIR=Y
    echo -n "Updating SCIFOR in $DIR_TARGET. Continue [Y/n]: "
    read _DIR;
    _DIR=`echo $_DIR |tr [:lower:] [:upper:]`
    [[ $_DIR == Y ]] || exit 1
else
    echo "Updating SCIFOR in $DIR_TARGET (quiet mode): "
    sleep 2
fi



cd $WRK_INSTALL


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
make update

if [ $? == 0 ];then
    make clean
    mv -vf $WRK_INSTALL/$LOG $DIR_TARGET/$LOG.update
    rm -vf $WRK_INSTALL/make.inc
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

cd $WRK_INSTALL
