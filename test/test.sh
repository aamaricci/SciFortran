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

while read DIR; do
    if [ -d $DIR ]; then
	if ls $DIR/*.x  1> /dev/null 2>&1; then
	    echo "TESTING $DIR"
	    cd $DIR
	    pwd
	    for exe in *.x
	    do
		echo "Running $exe:"
		./$exe
		echo ""
	    done
	    cd $HERE
	fi
    fi
done<list_dir
    

