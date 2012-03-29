#!/bin/bash

get_line()
{
if [ "$1" = "--help" -o "$1" = "-h" ];then
   echo "get_line: select given lines from file or pipe"	
   echo "cat ...|get_line 1 2 3... >"
   return
fi
ARG=
for LINE in "$@"
do
    ARG="$ARG${LINE}p;"
done
sed -n "$ARG"
}

pack()
{
if [ "$1" = "--help" -o "$1" = "-h" ];then
   echo "pack: remove blank lines from SINGLE number columns"	
   echo "cat ...|pack"
   echo "cat ...|awk '{print $X}'|pack"
   echo "pack < file"
   return
fi
sed '/^$/d'
}


cmplx()
{
if [ "$1" = "--help" -o "$1" = "-h" ];then
   echo "cmplx: format two columns in a fortran complex column:"	
   echo "cmplx [-e/--exchange] <col1> <col2>"
   return
fi
if [ "$1" = "-e" -o "$1" = "--exchange" ];then
    awk '{print "("$2","$1")"}'
else
    awk '{print "("$1","$2")"}'
fi
}

fcmplx()
{
if [ "$1" = "--help" -o "$1" = "-h" ];then
   echo "fcmplx: format a fortran complx in a real and imag parts:"	
   echo "fcmplx [-e/--exchange] <cmplx> "
   return
fi
if [ "$1" = "-e" -o "$1" = "--exchange" ];then
    sed -e 's/,/  /g' -e 's/(/ /g' -e 's/)/ /g'|awk '{print $2,$1}'
else
    sed -e 's/,/  /g' -e 's/(/ /g' -e 's/)/ /g'
fi
}

filter () 
{ 
if [ "$1" = "--help" -o "$1" = "-h" ];then
   echo "filter: use awk to filter number columns using input string"	
   echo "filter 'expr($1,$2,...)'"
   return
fi
EXPR=$1
awk '{print '"$EXPR"'}' 
}


reshape_histogram()
{
    if [ "$1" = "--help" -o "$1" = "-h" ];then
	echo "reshape_histogram: given histogram output, produce printable histogram "	
	echo "cat ...|reshape_histogram"
	echo "reshape_histogram < file(w/ histogram output)"
	return
    fi
    awk '{print $1, $3 ; print $2, $3}' 
}


get_histogram()
{
    if [ "$1" = "-h" -o "$1" = "--help" ];then
	echo "get_histogram: produce a printable histogram from a set of variables "	
	echo "cat ...|get_histogram"
	echo "get_histogram < file"
       return
    fi			
    histogram|awk '{print $1, $3 ; print $2, $3}' 
}



#FILE SIZE STORE LIMIT IN Kb (IOTOOLS)
export STORE_SIZE=2048
set_store_size()
{
    SIZE=$1
    export STORE_SIZE=$SIZE
    echo $STORE_SIZE
}


log_grid()
{
    if [ "$1" = "-h" -o "$1" = "--help" ];then
	echo "usage: $0 <1:log_base> <2:base_range> <3:num_steps>"
	echo "default 1:10"
	echo "default 2:1:3"
	echo "default 3:seq 0 9"
	return
    fi
    if [ -z $1 ];then
	BASE=10
    else
	BASE=$1
    fi
    if [ -z $2 ];then
	BRANGE=`seq 1 3`
    else
	BRANGE=`seq 1 $2`
    fi
    if [ -z $3 ];then
	NSTEP=`seq 0 9`
    else
	NSTEP=`seq 0 $3`
    fi
    for n in $BRANGE ; do 
     	DT=`echo "($BASE^$n-$BASE^($n-1))/$BASE" |bc -l`
     	for i in $NSTEP ; do 
     	    b=`echo "$BASE^($n-1)+$i*$DT"|bc -l` 
     	    echo $b 
     	done
    done
}
