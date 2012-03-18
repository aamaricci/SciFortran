#!/bin/bash

#INSTALLATION SCRIPT FOR THE SCIentific FORtran library.

#DEFINE SOME VARIABLES:
#==============================
SFROOT=`pwd`
SFDIR=$SFROOT/sf
SFLOCAL=$SFROOT/local
SFETC=$SFDIR/etc
SFSRC=$SFDIR/src
SFBIN=$SFDIR/bin
SFBLAS=$SFLOCAL/blas
SFLAPACK=$SFLOCAL/lapack
SFFFTW3=$SFLOCAL/fftw3




#START THE DIALOGUE WITH USERS:
#==============================
echo ""
echo "The SciFor library will be installed in dir $SFROOT "
echo "= = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ="
echo ""
echo ""
echo ""



#ASK WHICH COMPILER TO USE:
#==============================
echo "Please choose compiler? [gfortran(default),ifort] (others compiler not supported yet)"
read FCTMP
if [ -z $FCTMP ]; then
    FCTMP=gfortran
fi
if [ $FCTMP != "ifort" -o $FCTMP != "gfortran" ];then
    return 1
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
    echo "Here is a list of possible directories:"
    locate mklvars
    echo ""
    echo "Please write the absolute path to MKL tree directory"
    read MKLDIR
    echo ""
    echo "Please, write the aboslute path to the MKL vars script (e.g. mklvarsem64t.sh)"
    read MKLVARSH
    echo ""
    echo "MKL is installed in $MKLDIR"
    echo "The mkl vars script is $MKLVARSH"
else
    echo "MKL is not installed in your system"
fi
#add verification MKL is currently installed, perform a check
echo ""
echo ""
echo ""





#EDIT THE OPT.MK FILE:
#==============================
echo "Setup the opt.mk file, contains compilation options"
cd $SFETC
if [ $FCTMP == "gfortran" ];then
    cat <<EOF > opt.mk
OPT = -O3
STD = -O2 -fbounds-check
DEB = -O0 -g3 -fbounds-check -fbacktrace #-Wall -Wextra -Wconversion -pedantic
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
echo 'FFLAG += $(STD) -static' >> opt.mk
echo 'DFLAG += $(DEB) -static' >> opt.mk
cat opt.mk
sleep 2
echo ""
echo ""
echo ""

cd $SFROOT
#COMPILE THE BLAS/LAPACK/FFTW3 LIBRARIES
#=======================================
#start compiling the necessary libraries:
#blas:
echo "Compile local version of BLAS"
sleep 1
cd $SFBLAS
pwd
echo "output of the *make call is in make.log"
cp make.inc.$FCTMP make.inc
make 2>&1 > make.log 
DOBLAS=$?
echo "Success: $DOBLAS"
sleep 2
echo ""
echo ""
echo ""


# #lapack
echo "Compile local version of LAPACK"
sleep 1
cd $SFLAPACK
pwd
echo "output of the *make call is in make.log"
cp make.inc.$FCTMP make.inc
make 2>&1 > make.log
DOLAPACK=$?
echo "Success: $DOLAPACK"
sleep 2
echo ""
echo ""
echo ""

# #fftw3
echo "Compile local version of FFTW3"
cd $SFFFTW3
sleep 1
pwd
echo "output of the *./configure;make;make install call is in make.log"
./configure --prefix=`pwd` FC=$FCTMP 2>&1 > make.log
make 2>&1 >> make.log
make install 2>&1 >> make.log
DOFFTW3=$?
echo "Success: $DOFFTW3"
sleep 2
echo ""
echo ""
echo ""

cd $SFDIR
ln -sv $SFLOCAL $SFDIR/local



#EDIT THE LIBRARY.CONF FILE:
#==============================
echo "Setup the library.conf file"
cd $SFETC
cat <<EOF > library.conf
export FC=$FCTMP
export SFBLAS=$SFBLAS
export SFLAPACK=$SFLAPACK
export SFFFTW3=$SFFFTW3
EOF

if [ ! -z $MKLDIR ];then
    cat <<EOF >> library.conf
export MKLDIR=$MKLDIR
export MKLVARSH=$MKLVARSH
EOF
fi
cat library.conf
echo ""
echo ""
echo ""
sleep 2




#COMPILE THE LIBRARY
#=======================================
cd $SFSRC
if [ ! -z $MKLDIR ];then
    ln -sv FFT_MKL.f90 FFTGF.90
elif [ "$DOFFTW3" -eq "0" ];then
    ln -sv FFT_FFTW3.f90 FFTGF.90
else
    ln -sv FFT_NR.f90 FFTGF.90
fi

source $SFDIR/bin/mylibvars.sh
make 
echo""
echo""
echo""


echo 'Please be sure to add the following line to your shell init file (e.g. .bashrc, .profile, .bash_profile,etc..)'
echo "export SFDIR=$SFDIR"
echo 'source $SFDIR/bin/mylibvars.sh'

echo ""
echo "done. good bye"
echo "Please open a new terminal session"
echo ""
echo ""


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









