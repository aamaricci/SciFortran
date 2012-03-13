#!/bin/bash
HERE=`pwd`
LIST=`ls -l |egrep '^d' |awk '{print $8}'`

# for DIR in $LIST
# do
#     cd $HERE/$DIR
#     cat Makefile |sed '/include/d' > tmp
#     echo 'include $(SFDIR)/etc/opt.mk' > Makefile
#     cat tmp >> Makefile
#     rm -fv tmp
# done
# exit

EXCLUDED='COMVARS FFT_MKL FFT_NR FFT_FFTW3 FFT_GSL' 



for DIR in $EXCLUDED
do 
    LIST=`echo $LIST | sed "s/\b$DIR\b//g"`
done

LIST='COMVARS '`echo $LIST`' FFT'
echo $LIST
cd $HERE
for DIR in $LIST
do 
    echo "===========$DIR============="
    cd $HERE/$DIR
    make  ; FAIL=$?
    if [ "$FAIL" -ne "0" ];then
	echo "####################"
	echo "ERROR UPDATING $DIR!"
	echo "####################"
	ERRLIST=`echo $ERRLIST`" $DIR"
    fi
    echo ''
done
cd $HERE

if [ "$ERRLIST" = "" ];then 
    echo "all libraries have been updated"
else
    echo "the following libraries could not be updated:"
    echo $ERRLIST
fi
