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

cp -v src/inputAHM.conf bin/
cd bin/

for exe in *.x
do
    ./$exe 
done
    
    

