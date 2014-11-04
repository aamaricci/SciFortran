#!/bin/bash

#THIS ADD A GIVEN LIBRARY TO THE SYSTEM:
add_library_to_system(){
    LIB=$1
    if [ ! -d $LIB ];then
	echo "Error: can not find $libdir directory" >&2
	exit
    else
	echo "Add Scientific Fortran (SciFor) library" >&2
	#
	if [ -d $LIB/lib ];then
	    if [ -z "${LD_LIBRARY_PATH}" ];then
		export LD_LIBRARY_PATH=$LIB/lib
	    else
		export LD_LIBRARY_PATH=$LIB/lib:$LD_LIBRARY_PATH
	fi
	#
	    if [ -z "${LIBRARY_PATH}" ];then
		export LIBRARY_PATH=$LIB/lib
	    else
		export LIBRARY_PATH=$LIB/lib:$LIBRARY_PATH
	    fi
	fi
	#
	if [ -d $LIB/include ];then
	    if [ -z "${FPATH}" ];then
		export FPATH=$LIB/include
	    else
		export FPATH=$LIB/include:$FPATH
	    fi
	#
	    if [ -z "${CPATH}" ];then
		export CPATH=$LIB/include
	    else
		export CPATH=$LIB/include:$CPATH
	    fi
	fi
	#
	if [ -d $LIB/man ];then
	    if [ -z "${MANPATH}" ];then
		export MANPATH=$LIB/man
	    else
		export MANPATH=$LIB/man:$MANPATH
	    fi
	fi
	#
	if [ -d $LIB/bin ];then
	    if [ -z "${PATH}" ];then
		export PATH=$LIB/bin
	    else
		export PATH=$LIB/bin:$PATH
	    fi
	fi
    fi
}


#SETUP SCIFOR LIBRARY:
if [ -z $1 ] || [ -z $2 ]; then
    echo "Usage: $0 <path_to_library_root> <platform [intel,gnu,intel_debug,gnu_debug>" >&2
    return 1
fi

root=$1
plat=$2
libdir=$root/$plat

if [ ! -d $root ];then
    echo "Error: can not find $root directory" >&2
fi



add_library_to_system $libdir
export SFROOT=$libdir
