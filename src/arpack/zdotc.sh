#!/bin/bash
ARGUMENTS=${@}
PATH_TO_ME=${0}
SCRIPT=${0##*/}
LOGFILE=log.${SCRIPT%.sh}
>$LOGFILE
exec >  >(tee -a $LOGFILE)
exec 2> >(tee -a $LOGFILE >&2)

PATCH=zdotc.patch
#TEST PATCH:
printf "Testing $PATCH..."
patch -Rsfp0 --dry-run < $PATCH
printf "Done\n"
if [ $? -eq 0 ];then
    echo "Patch applied"
else
    echo "Patch not applied... applying $PATCH"
    patch -s -p0 < $PATCH
    if [ $? -eq 0 ];then
	echo "Patch successfully applied"
	return 0
    else
	echo "ERROR in applying $PATCH"
	return 1
    fi
fi
