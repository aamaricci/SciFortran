#!/bin/bash
NAME=SCIFOR
LIBNAME=libscifor.a

################### >>>> FUNCTIONS <<<<<< ###################
usage(){
    echo "usage:  $0 (-h,--help) [plat] (opt) "
    echo "   -h,--help: this help"
    echo "        plat: (variable) + intel,intel_debug,gnu,gnu_debug - specifies the actual platform/compiler to use"
    echo "         opt: (variable) + arpack blas fftpack lapack minpack quadpack dmftt - to install a single 3rd party lib. "
    echo "                         + clean - to clean out the former compilation for the specified plat." 
    echo "                         + wdmftt - to complete SciFor with DMFT_TOOLS library."
    exit
}



if [ -z $1 ] || [ $1 == "-h" ] || [ $1 == "--help"  ];then
    usage
fi

LIST_FC="gnu gnu_debug intel intel_debug"
[[ $LIST_FC =~ (^|[[:space:]])"$1"($|[[:space:]]) ]] || usage

LIST_OPT="clean arpack blas fftpack lapack minpack quadpack dmftt wdmftt"
OPT=""
if [ ! -z $2 ];then
    OPT=$2
    [[ "$LIST_OPT" =~ (^|[[:space:]])"$OPT"($|[[:space:]]) ]] || usage
fi

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
OPT_INSTALL=$WRK_INSTALL/opt
ENVMOD_INSTALL=$ETC_INSTALL/environment_modules
SRC_INSTALL=$WRK_INSTALL/src
PREFIX=$(nparent_dir $WRK_INSTALL 1)

OPT_LIBS="arpack blas fftpack lapack minpack quadpack"

DIR_TARGET=$PREFIX/$PLAT
BIN_TARGET=$DIR_TARGET/bin
ETC_TARGET=$DIR_TARGET/etc
LIB_TARGET=$DIR_TARGET/lib
INC_TARGET=$DIR_TARGET/include

#test MPI
FMPI=mpif90
which $FMPI >/dev/null 2>&1
WMPI=$?

create_makeinc(){    
    local PLAT=$1
    cd $WRK_INSTALL
    echo "Creating directories:" >&2
    mkdir -pv $DIR_TARGET
    mkdir -pv $BIN_TARGET
    mkdir -pv $ETC_TARGET/modules/$LNAME
    mkdir -pv $LIB_TARGET
    mkdir -pv $INC_TARGET
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

    echo "Copying init script for $UNAME" >&2
    cp -fv $BIN_INSTALL/scifor_completion.sh  $BIN_TARGET/scifor_completion.sh
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
    cat $ENVMOD_INSTALL/module >> $ETC_TARGET/modules/${LNAME}/$PLAT
    echo "" >&2
    echo "Compiling $UNAME library on platform $PLAT:">&2
    echo "" >&2
}


############# >>>> WORK STARTS HERE <<<<<< #################
HERE=$(pwd)

#CLEAN ALL IF YOU ARE REQUIRED TO:
case $OPT in
    arpack)
	if [ $WMPI == 1 ];then
	    cd $OPT_INSTALL/arpack
	    ./install.sh $PLAT
	    exit
	else
	    cd $OPT_INSTALL/parpack
	    ./install.sh $PLAT
	    exit
	fi
	;;
    blas)
	cd $OPT_INSTALL/blas
	./install.sh $PLAT
	exit
	;;
    dmftt)
	cd $OPT_INSTALL/dmft_tools
	./install.sh $PLAT
	exit
	;;
    fftpack)
	cd $OPT_INSTALL/fftpack
	./install.sh $PLAT
	exit
	;;
    lapack)
	cd $OPT_INSTALL/lapack
	./install.sh $PLAT
	exit
	;;
    minpack)
	cd $OPT_INSTALL/minpack
	./install.sh $PLAT
	exit
	;;
    quadpack)
	cd $OPT_INSTALL/quadpack
	./install.sh $PLAT
	exit
	;;
    clean)
	echo "Cleaning the .o files in all OPT libraries for platform $PLAT: " >&2
	cd $OPT_INSTALL
	for dir in *;do
	    cd $dir/
	    ./install.sh $PLAT clean
	    cd ../
	done
	exit
	;;
esac

# if [ $OPT == "clean" ];then
#     echo "Cleaning the .o files in all OPT libraries for platform $PLAT: " >&2
#     cd $OPT_INSTALL
#     for dir in *;do
# 	cd $dir/
# 	./install.sh $PLAT clean
# 	cd ../
#     done
#     exit
# fi

cd $HERE

#INSTALLING OPT LIBRARIES IF NOT PRESENT IN THE TARGET DIRECTORY:
for LIB in $OPT_LIBS;do
    case $LIB in
	arpack)
	    if [ $WMPI == 1 ];then
		if [ ! -e $LIB_TARGET/libarpack.a ];then
		    cd $OPT_INSTALL/arpack
		    ./install.sh $PLAT
		else
		    echo "$LIB_TARGET/libarpack.a exists. skip" >&2
		    continue
		fi
	    else
		if [ ! -e $LIB_TARGET/libparpack.a ] ;then
		    cd $OPT_INSTALL/parpack
		    ./install.sh $PLAT
		else
		    echo "$LIB_TARGET/libarpack.a exists. skip" >&2
		    echo "$LIB_TARGET/libparpack.a exists. skip" >&2
		    continue
		fi
	    fi
            ;;
	blas)
	    if [ ! -e $LIB_TARGET/libblas.a ];then
		cd $OPT_INSTALL/$LIB
		./install.sh $PLAT
	    else
		echo "$LIB_TARGET/libblas.a exists. skip" >&2
	    fi
	    ;;
	fftpack)
	    if [ ! -e $LIB_TARGET/libfftpack.a ];then
		cd $OPT_INSTALL/$LIB
		./install.sh $PLAT
	    else
		echo "$LIB_TARGET/libfftpack.a exists. skip" >&2
	    fi
	    ;;
	lapack)
	    if [ ! -e $LIB_TARGET/liblapack.a ];then
		cd $OPT_INSTALL/$LIB
		./install.sh $PLAT
	    else
		echo "$LIB_TARGET/liblapack.a exists. skip" >&2
	    fi
	    ;;
	minpack)
	    if [ ! -e $LIB_TARGET/libminpack.a ];then
		cd $OPT_INSTALL/$LIB
		./install.sh $PLAT
	    else
		echo "$LIB_TARGET/libminpack.a exists. skip" >&2
	    fi
	    ;;
	quadpack)
	    if [ ! -e $LIB_TARGET/libquadpack.a ];then
		cd $OPT_INSTALL/$LIB
		./install.sh $PLAT
	    else
		echo "$LIB_TARGET/libquadpack.a exists. skip" >&2
	    fi
	    ;;
    esac	
done
cd $HERE


#CREATING THE SCIFOR LIBRARY:
rm -fv $LIB_TARGET/$LIBNAME
create_makeinc $PLAT
make all
if [ $? == 0 ];then
    make clean
    mv -vf $WRK_INSTALL/make.inc $ETC_TARGET/make.inc.scifor
else
    echo "Error from Makefile. STOP here."
    exit 1
fi


#IF REQUIRED GENERATING THE DMFT_TOOLS LIBRARY FROM OPT
if [ ! -z $2 ] && [ "$OPT" == "wdmftt" ];then
    HERE=$(pwd)
    echo "Generating the DMFT_TOOLS library for $PLAT" >&2
    cd $OPT_INSTALL/dmft_tools
    ./install.sh $PLAT
    cd $HERE
fi


#LAST TOUCH COPY THE CONFIGVARS AND CREATE THE USER MODULES FILE. 
#PRINT USAGE DETAILS.
CONFIGFILE=$PREFIX/$PLAT/bin/configvars.sh
MODULEFILE=$PREFIX/$PLAT/etc/modules/$LNAME/$PLAT
mkdir -pv $HOME/.modules.d/$LNAME
cp -vf $MODULEFILE $HOME/.modules.d/$LNAME/$PLAT
echo "" >&2
echo "USAGE:" >&2
echo "" >&2
echo "  >  source $CONFIGFILE" >&2
echo "(place such line into your bash profile [e.g. .bashrc])" >&2
echo "">&2
module avail >/dev/null 2>&1 
if [ $? == 0 ];then
    echo "OR load the $UNAME modules: (place the following lines into your bash profile [e.g. .bashrc])" >&2
    echo "  > module use $HOME/.modules.d/applications" >&2
    echo "  > module load $LNAME/$PLAT" >&2
    echo "(place these lines into your bash profile [e.g. .bashrc])" >&2
    echo "">&2
fi
echo "">&2
echo "Enjoy... (for additional info: adriano.amaricciATgmail.com)">&2
echo "">&2
