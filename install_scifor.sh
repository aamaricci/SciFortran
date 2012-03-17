#!/bin/bash

#INSTALLATION SCRIPT FOR THE SCIentific FORtran library.

#DEFINE SOME VARIABLES:
#==============================
HERE=`pwd`
export SFDIR=$HERE/scifor
SFLOCAL=$HERE/local
SFETC=$SFDIR/etc
SFSRC=$SFDIR/src
SFBIN=$SFDIR/bin
SFBLAS=$SFLOCAL/blas
SFLAPACK=$SFLOCAL/lapack
SFFFTW3=$SFLOCAL/fftw3




#START THE DIALOGUE WITH USERS:
#==============================
echo ""
echo "The SciFor library will be installed in dir $SFDIR "
echo "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="
echo ""
echo ""
echo ""




#ASK WHICH COMPILER TO USE:
#==============================
echo "Please choose compiler? [gfortran*,ifort,xlf,...]"
read FCTMP
if [ -z $FCTMP ]; then
    FCTMP=gfortran
fi
echo "You choose to use **$FCTMP** compiler"
echo ""
echo ""
echo ""




#CHECK PRESENCE OF MKL
#==============================
echo "Do you have MKL installed in your system? [N,y]"
read MKLANSWER
if [ -z $MKLANSWER ];then
    MKLANSWER=n
fi
if [ $MKLANSWER == "y" -o $MKLANSWER == "Y" ];then
    echo "Please write the absolute path to MKL tree directory"
    read MKLDIR
    echo "Please, write the aboslute path to the MKL vars script (e.g. mklvarsem64t.sh)"
    read MKLVARSH
    echo "MKL is installed in $MKLDIR"
    echo "The mkl vars script is $MKLVARSH"
else
    echo "MKL is not installed in your system"
fi
#add verification MKL is currently installed, perform a check
echo ""
echo ""
echo ""







#EDIT THE LIBRARY.CONF FILE:
#==============================
echo "Setup the library.conf file"
cd $SFETC
cat <<EOF > library.conf
export FC=$FCTMP
export FFTW3DIR=$SFDIR/local/fftw3
EOF

if [ ! -z $MKLDIR ];then
    cat <<EOF >> library.conf
export MKLDIR=$MKLDIR
source $MKLVARSH
EOF
fi
echo ""
echo ""
echo ""





#EDIT THE OPT.MK FILE:
#==============================
echo "Setup the opt.mk file, contains compilation options"
if [ $FCTMP == "gfortran" ];then
    cat <<EOF > opt.mk
OPT = -pg -O3
STD = -pg -O1
DEB = -pg -O0 -g3 -fbounds-check -fbacktrace #-Wall -Wextra -Wconversion -pedantic
EOF
elif [ $FCTMP == "ifort" ];then
    cat <<EOF > opt.mk
OPT =  -O3 -ftz -assume nobuffered_io -openmp #-parallel
STD =  -O2 -assume nobuffered_io
DEB =  -p -traceback -O0 -g -debug -fpe0 -traceback  #-static-intel -check all
EOF
else
    cat <<EOF >opt.mk
OPT=
STD=
DEB=
EOF
fi


cat <<EOF >> opt.mk
#########################################################################
FFLAG += $(STD) -static
DFLAG += $(DEB) -static
EOF
echo ""
echo ""
echo ""
cd $SFDIR


#COMPILE THE BLAS/LAPACK/FFTW3 LIBRARIES
#=======================================
#start compiling the necessary libraries:
#blas:
echo "Compile local version of BLAS"
cd $SFBLAS
##
pwd
##
# cp make.inc.$FCTMP make.inc
# make 2>&1 > make.log 


# #lapack
echo "Compile local version of LAPACK"
cd $SFLAPACK
##
pwd
##
# cp make.inc.$FCTMP make.inc
# make 2>&1 > make.log


# #fftw3
echo "Compile local version of FFTW3"
cd $SFFFTW3
##
pwd
##
# ./configure --prefix=`pwd` 2>&1 > make.log
# make 2>&1 >> make.log
# make install 2>&1 >> make.log
#GET THE BOOLEAN $? ouput: check if compilation was ok or not, needed to setup FFT
EFFTW3=1
# cd $SFDIR
# echo ""
# echo ""
# echo ""

cd $SFDIR
ln -s $SFLOCAL $SFDIR/local

#COMPILE THE LIBRARY
#=======================================
cd $SFSRC
if [ ! -z $MKLDIR ];then
    rm -fr FFT
    ln -s ./FFT_MKL ./FFT
elif [ $EFFTW3 == "0" ];then
    rm -fr FFT
    ln -s ./FFT_FFTW3 ./FFT
else
    rm -fr FFT
    ln -s ./FFT_NR ./FFT
fi

sh update_lib.sh





##################################################################
###REPO
# #CHECK PRESENCE OF GSL/FGSL
# #==============================
# echo "Do you have GSL/FGSL installed in your system? [N,y]"
# read GSLANSWER
# if [ -z $GSLANSWER ];then
#     GSLANSWER=n
# fi
# if [ $GSLANSWER == "y" -o $GSLANSWER == "Y" ];then
#     echo "Please write the absolute path to GSL tree directory"
#     read GSLDIR
#     echo "Please write the absolute path to FGSL tree directory"
#     read FGSLDIR
#     echo "GSL is installed in $GSLDIR"
#     echo "FGSL is installed in $FGSLDIR"
# else
#     echo "GSL/FGSL is not installed in your system"
# fi
# #add verification GSL/FGSL is currently installed, perform a check
# echo ""
# echo ""
# echo ""

# if [ ! -z $GSLDIR ];then
#     cat <<EOF >> library.conf
# export GSLDIR=$GSLDIR
# EOF
# fi

# if [ ! -z $FGSLDIR ];then
#     cat <<EOF >> library.conf
# export FGSLDIR=$FGSLDIR
# EOF
# fi









