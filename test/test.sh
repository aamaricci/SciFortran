#Code to run all tests in test/bin/
#N.B. test should end by .x

set -e

if [ ! -d "bin" ]
then
    echo "\e[31m ERROR \e[0m"
    echo " There is no *bin* directory"
    echo " Try  'make all' before testing"
    return 1
fi

cd bin/


cat <<EOF > SF_PARSE_INPUT/input.conf
DBLE=0.5000000000                           !comment A dble
L=1000                                      !An integer
BOOL=F                                      !A boolean
DBLE_VEC=1d0,1d0,0d0			    !A dble vector
STRING=input.conf			    !A string 
EOF

HERE=`pwd`
echo $HERE
for dir in {SF_DERIVATE,SF_FFT,SF_INIT,SF_INTEGRATE,SF_INTERPOLATE,SF_IOTOOLS,SF_LINALG,SF_MPI,SF_OPTIMIZE,SF_PARSE_INPUT,SF_RANDOM,SF_SPECIAL,SF_SP_LINALG};
do
    if compgen -G "$dir/*.x" > /dev/null ;then
	echo "TESTING $dir"
	cd $dir
	pwd
	for exe in *.x
	do
	    echo "Running $exe:"
	    ./$exe
	    echo ""
	done
	cd $HERE
    fi
done
    

