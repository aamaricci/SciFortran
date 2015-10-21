#!/bin/bash
if [ -z "$1" ];then
    echo "target executable not defined. stop."
    exit 1
fi
BRANCH=_$(git rev-parse --abbrev-ref HEAD)
if [ $BRANCH == "_master" ];then
    BRANCH=
fi
SUFFIX=$BRANCH
if [ -z "$BRANCH" ];then
    SUFFIX=
fi
DIR_BASH_COMPLETION=$HOME/.bash_completion.d
mkdir -p $DIR_BASH_COMPLETION
EXE_DRIVER=$1
OPT_INPUT=`grep -e parse_input -e parse_cmd *.f90 $EXE_DRIVER | cut -d "(" -f2 | cut -d ")" -f1|cut -d"," -f1|awk '{print $1"="}'|tr "\\n" ' '`
EXE=`basename $EXE_DRIVER`
EXE=${EXE%.f90}
cat <<EOF > $DIR_BASH_COMPLETION/$EXE
_${EXE}${SUFFIX}()
{
local cur=\${COMP_WORDS[COMP_CWORD]}
COMPREPLY=( \$(compgen -W "${OPT_INPUT}" \${cur}) )
}
complete -F _${EXE}${SUFFIX} ${EXE}${SUFFIX}
EOF
chmod u+x $DIR_BASH_COMPLETION/$EXE
source $DIR_BASH_COMPLETION/$EXE

for PROFILE in .bash_profile .bashrc .profile
do
    if [ -e $PROFILE ];then
	sed  -n -e  '/#BEGIN_SCIFOR/{p; :a; N; /#END_SCIFOR/!ba; s/.*\n//}; p' $PROFILE |sed '/#BEGIN_SCIFOR/d' |sed '/#END_SCIFOR/d' > tmp
	cp -vf $PROFILE $PROFILE.orig
	mv -vf tmp $PROFILE
	cat <<EOF >> $HOME/$PROFILE
#BEGIN_SCIFOR COMPLETION
if [ -d $DIR_BASH_COMPLETION ];then
   for file in $DIR_BASH_COMPLETION/*
   do
    . $file
   done
# FIXME: 
#END_SCIFOR COMPLETION
EOF
    fi
done
