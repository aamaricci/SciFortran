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
    echo "$0  -p,--plat=FC_PLAT --prefix=PREFIX_DIR  [-o,--opt-lib=OPT_LIB --mpi=MPI_ROOT  -q,--quiet -c,--clean  -d,--debug  -h,--help]"
    echo ""
    echo "mandatory arguments:" 
    echo "    -p,--plat   : specifies the actual platform/compiler to use [default: intel] (valid: intel, gnu)"
    echo "    --prefix    : specifies the target directory [default: /opt/scifor]"
    echo ""    
    echo "optional arguments:" 
    echo "    -o,--opt-lib: to install a single 3rd party lib [arpack blas fftpack lapack minpack quadpack]. "
    echo "    --mpi       : specifies the MPI ROOT directory [default is the detected one]"
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
params="$(getopt -n "$0" --options p:o:qcdh --longoptions plat:,prefix:,mpi:,opt-lib:,quiet,clean,debug,help -- "$@")"
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
DEBUG=1
CLEAN=1
IMPI=1
OPT=""
QUIET=0
IDEFAULT=0
#DEFAULT PREFIX
PREFIX=/opt/scifor
#DEFAULT PLATFORM
PLAT=intel
VERSION=$(git describe --tags 2>/dev/null)
WRK_INSTALL=$(pwd)
BIN_INSTALL=$WRK_INSTALL/bin
ETC_INSTALL=$WRK_INSTALL/etc
OPT_INSTALL=$WRK_INSTALL/opt
ENVMOD_INSTALL=$ETC_INSTALL/environment_modules
SRC_INSTALL=$WRK_INSTALL/src


#>>> A COUPLE OF LISTS OF ALLOWED OPTIONS
LIST_FC="gnu intel"
LIST_OPT="parpack blas fftpack lapack minpack quadpack"


#>>> GO THROUGH THE INPUT ARGUMENTS. FOR EACH ONE IF REQUIRED TAKE ACTION BY SETTING VARIABLES.
while true
do
    case $1 in
	-p|--plat)
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
	-o|--opt-lib)
	    OPT=$2;
	    shift 2
	    [[ ! "$LIST_OPT" =~ (^|[[:space:]])"$OPT"($|[[:space:]]) ]] && {
		echo "Incorrect OPTional Library: $OPT";
		echo " available choices are: $LIST_OPT";
		exit 1
	    }
	    ;;
	--mpi)
	    IMPI=1;
	    MPIROOT=$2;
	    shift 2
	    ;;
	-q|--quiet) QUIET=1;shift ;;
	-c|--clean) CLEAN=0;shift ;;
	-d|--debug) DEBUG=0;shift ;;
        -h|--help) usage ;;
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

#>>> SET STANDARD NAMES FOR THE TARGET DIRECTORY
DIR_TARGET=$PREFIX/$PLAT
BIN_TARGET=$DIR_TARGET/bin
ETC_TARGET=$DIR_TARGET/etc
LIB_TARGET=$DIR_TARGET/lib
INC_TARGET=$DIR_TARGET/include
DIR_TARGET_W=0
if [ ! -w $PREFIX ];then
    DIR_TARGET_W=1
    echo "Can not create $DIR_TARGET: $PREFIX has no write access"
    sleep 0.5
    echo "Try to grant root privileges to create $DIR_TARGET for $USER:$GROUP"
    sudo -v
fi



#>>> TEST MPI
FMPI=mpif90
which $FMPI >/dev/null 2>&1
WMPI=$?
if [ $WMPI == 0 ];then
    _MPI=Y
    echo "Using the following MPI compiler:"
    mpif90 -show
    sleep 0.5
    echo "On platform: $PLAT"
    sleep 0.5
    MPIFC=$(mpif90 -show | awk '{print $1}')
    if [ "$MPIFC" != "$_FC_" ];then
	echo "Possible mismatch between MPI:+$MPIFC and Platform $PLAT:+$_FC_" 
	sleep 2
	echo "Fall back into the non-MPI mode"
	sleep 0.5
	WMPI=1
	# QUIET=0
    fi
    # if [ $QUIET == 0 ];then
    # 	echo -n "Continue (this will install P-Arpack) [Y/n]: "
    # 	read _MPI;
    # 	sleep 0.5
    # 	_MPI=`echo $_MPI |tr [:lower:] [:upper:]`
    # 	[[ $_MPI == Y ]] || exit 1
    # fi
fi



#>>> CREATE THE MAKE.INC IN THE FOLLOWING FUNCTION
create_makeinc(){
    local PLAT=$1
    cd $WRK_INSTALL
    case $PLAT in
	intel)
	    FC=ifort
	    FFLAGS="-O2 -static-intel"
	    MOPT="-module "
	    ;;
	gnu)
	    FC=gfortran
	    FFLAGS="-O2 -static"
	    MOPT=-J
	    ;;
	intel_debug)
	    FC=ifort
	    FFLAGS="-p -O0 -g -debug -fpe0 -traceback -check all,noarg_temp_created -static-intel"
	    MOPT="-module "
	    ;;
	gnu_debug)
	    FC=gfortran
	    FFLAGS="-O0 -p -g -Wall -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -static"
	    MOPT=-J
	    ;;
	ibm)
	    FC=xlf90
	    FFLAGS="-O2 -qarch=qp -qtune=qp"
	    MOPT="-qmoddir="
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
INC_TARGET=$INC_TARGET
LIB_SCIFOR=$LIB_TARGET/$LIBNAME
EOF
}


create_makeinc $PLAT
sleep 0.5
if [ $CLEAN == 0 ];then
    echo "Cleaning the .o files in all OPT libraries for platform $PLAT: " 
    cd $OPT_INSTALL
    for dir in *;do
	cd $dir/
	./install.sh $LIST_ARGS --clean
	cd ../
    done
    cd $HERE
    make clean
    exit 0
fi



if [ $QUIET == 0 ];then
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
sleep 0.5
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
sleep 0.5




#INSTALL A SINGLE OPTIONAL LIB IF REQUIRED:
case $OPT in
    arpack)
	if [ $WMPI == 1 ];then
	    cd $OPT_INSTALL/parpack
	    ./install.sh $LIST_ARGS
	    exit
	else
	    cd $OPT_INSTALL/parpack
	    ./install.sh $LIST_ARGS --mpi
	    exit
	fi
	;;
    blas)
	cd $OPT_INSTALL/blas
	./install.sh $LIST_ARGS
	exit
	;;
    fftpack)
	cd $OPT_INSTALL/fftpack
	./install.sh $LIST_ARGS
	exit
	;;
    lapack)
	cd $OPT_INSTALL/lapack
	./install.sh $LIST_ARGS
	exit
	;;
    minpack)
	cd $OPT_INSTALL/minpack
	./install.sh $LIST_ARGS
	exit
	;;
    quadpack)
	cd $OPT_INSTALL/quadpack
	./install.sh $LIST_ARGS
	exit
	;;
esac
cd $HERE


#CHECK THAT EACH OPT LIBRARY IS PRESENT IN THE INSTALLATION. IF NOT INSTALL IT
for LIB in $LIST_OPT;do
    case $LIB in
	arpack)
	    if [ $WMPI == 0 ];then
		if [ -e $LIB_TARGET/libarpack.a ] && [ -e $LIB_TARGET/libparpack.a ];then
		    echo "$LIB_TARGET/libarpack.a exists. skip" 
		    echo "$LIB_TARGET/libparpack.a exists. skip" 
		else
		    cd $OPT_INSTALL/parpack
		    ./install.sh $LIST_ARGS --mpi
		    [[ $? == 0 ]] || exit 1
		fi
	    else
		if [ -e $LIB_TARGET/libarpack.a ] ;then
		    echo "$LIB_TARGET/libarpack.a exists. skip" 
		else
		    cd $OPT_INSTALL/parpack
		    ./install.sh $LIST_ARGS
		    [[ $? == 0 ]] || exit 1
		fi
	    fi
            ;;
	blas)
	    if [ -e $LIB_TARGET/libblas.a ];then
		echo "$LIB_TARGET/libblas.a exists. skip" 
	    else
		cd $OPT_INSTALL/$LIB
		./install.sh $LIST_ARGS
		[[ $? == 0 ]] || exit 1
	    fi
	    ;;
	fftpack)
	    if [ -e $LIB_TARGET/libfftpack.a ];then
		echo "$LIB_TARGET/libfftpack.a exists. skip" 
	    else
		cd $OPT_INSTALL/$LIB
		./install.sh $LIST_ARGS
		[[ $? == 0 ]] || exit 1
	    fi
	    ;;
	lapack)
	    if [ -e $LIB_TARGET/liblapack.a ];then
		echo "$LIB_TARGET/liblapack.a exists. skip" 
	    else
		cd $OPT_INSTALL/$LIB
		./install.sh $LIST_ARGS
		[[ $? == 0 ]] || exit 1
	    fi
	    ;;
	minpack)
	    if [ -e $LIB_TARGET/libminpack.a ];then
		echo "$LIB_TARGET/libminpack.a exists. skip" 
	    else
		cd $OPT_INSTALL/$LIB
		./install.sh $LIST_ARGS
		[[ $? == 0 ]] || exit 1
	    fi
	    ;;
	quadpack)
	    if [ -e $LIB_TARGET/libquadpack.a ];then
		echo "$LIB_TARGET/libquadpack.a exists. skip" 
	    else
		cd $OPT_INSTALL/$LIB
		./install.sh $LIST_ARGS
		[[ $? == 0 ]] || exit 1
	    fi
	    ;;
    esac	
done
cd $HERE


echo "Copying init script for $UNAME" 
cp -fv $BIN_INSTALL/scifor_completion.sh  $BIN_TARGET/scifor_completion.sh
cp -fv $BIN_INSTALL/configvars.sh $BIN_TARGET/configvars.sh
cat <<EOF >> $BIN_TARGET/configvars.sh
add_library_to_system ${PREFIX}/${PLAT}
EOF
echo "" 
sleep 0.5


echo "Generating environment module file for $UNAME" 
cat <<EOF > $ETC_TARGET/modules/$LNAME/$PLAT
#%Modules
set	root	$PREFIX
set	plat	$PLAT
set	version	"$VERSION ($PLAT)"
EOF
cat $ENVMOD_INSTALL/module >> $ETC_TARGET/modules/${LNAME}/$PLAT
echo "" 
sleep 0.5

echo "Compiling $UNAME library on platform $PLAT:"
echo "" 
sleep 0.5



#CREATING THE SCIFOR LIBRARY:
rm -fv $LIB_TARGET/$LIBNAME
make all
if [ $? == 0 ];then
    make clean
    mv -vf $WRK_INSTALL/make.inc $ETC_TARGET/make.inc.scifor
else
    echo "Error from Makefile. STOP here."
    exit 1
fi
sleep 0.5


#LAST TOUCH COPY THE CONFIGVARS AND CREATE THE USER MODULES FILE. PRINT USAGE DETAILS.
CONFIGFILE=$PREFIX/$PLAT/bin/configvars.sh
MODULEFILE=$PREFIX/$PLAT/etc/modules/$LNAME/$PLAT
mkdir -pv $HOME/.modules.d/$LNAME
cp -vf $MODULEFILE $HOME/.modules.d/$LNAME/$PLAT

echo ""
echo "Merge all libraries into ${LIBNAME%.a}_all.a"
sleep 0.5
cd $LIB_TARGET
rm -fv ${LIBNAME%.a}_all.a
for LIB in lib*.a;
do 
    DIRLIB=tmp_${LIB%.a};
    mkdir -pv $DIRLIB;
    ar -xv $LIB && mv -v *.o $DIRLIB/;
done
ar -crv ${LIBNAME%.a}_all.a `ls tmp_*/*.o |sort |uniq` && rm -rfv tmp_*
sleep 0.5
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
echo "   $ module use $HOME/.modules.d/applications" 
echo "   $ module load $LNAME/$PLAT" 
echo "(or copy these lines into your bash profile [e.g. .bashrc])" 
echo ""
fi
echo ""
echo "Enjoy... (for info: adriano.amaricciATgmail.com)"
echo ""
cd $HERE
cp -v $LOG $PREFIX/$PLAT/
