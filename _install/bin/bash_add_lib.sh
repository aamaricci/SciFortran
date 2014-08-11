#!/bin/bash

add_library_to_system(){
    LIB=$1
    if [ ! -d $LIB ];then
	echo "directory $LIB does not exist: skip"
    else
	echo "Adding $LIB to system..."
	#
	if [ -z "${LD_LIBRARY_PATH}" ];then
	    export LD_LIBRARY_PATH=${LIB}/lib
	else
	    export LD_LIBRARY_PATH=${LIB}/lib:$LD_LIBRARY_PATH
	fi
	#
	if [ -z "${LIBRARY_PATH}" ];then
	    export LIBRARY_PATH=${LIB}/lib
	else
	    export LIBRARY_PATH=${LIB}/lib:$LIBRARY_PATH
	fi
	#
	if [ -z "${FPATH}" ];then
	    export FPATH=${LIB}/include
	else
	    export FPATH=${LIB}/include:$FPATH
	fi
	#
	if [ -z "${CPATH}" ];then
	    export CPATH=${LIB}/include
	else
	    export CPATH=${LIB}/include:$CPATH
	fi
	#
	if [ -z "${MANPATH}" ];then
	    export MANPATH=${LIB}/man
	else
	    export MANPATH=${LIB}/man:$MANPATH
	fi
    fi
}

if [ -z $1 ];then
    usage "$0 <path to library to be installed>"
    exit
fi
add_library_to_system $1
