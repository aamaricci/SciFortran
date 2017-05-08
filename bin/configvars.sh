#!/bin/bash
#THIS ADD A GIVEN LIBRARY TO THE SYSTEM:
add_library_to_system(){
    LIB=$1
    if [ ! -d $LIB ];then
	echo "Error: can not find $LIB directory" >&2
    else
	echo "Add $LIB to env" >&2
	#
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
	#
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
	if [ -z "${MANPATH}" ];then
	    export MANPATH=$LIB/man
	else
	    export MANPATH=$LIB/man:$MANPATH
	fi
	if [ -z "${PATH}" ];then
	    export PATH=$LIB/bin
	else
	    export PATH=$LIB/bin:$PATH
	fi
	if [ -z "${PKG_CONFIG_PATH}" ];then
	    export PKG_CONFIG_PATH=$HOME/.pkgconfig
	else
	    export PKG_CONFIG_PATH=$HOME/.pkgconfig:$PKG_CONFIG_PATH
	fi
    fi
}
