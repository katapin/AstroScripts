#!/bin/sh
#
#Wrapper for xmm-sas programs called from scripts
#It executes "sasinit" before run the program
#Synax: xmmcall.sh string_for_call


if [ "x${1}" == "x" ]
then
    echo "Error: nothing to call"
    exit 1;
fi

PROGNAME="$1"       #Program name
shift
ARGS=("$@")         #Array of arguments

#shopt -s expand_aliases
source ~/.bashrc
sasinit >/dev/null

Cwhite="\e[1;37m"
CNC="\e[0m";
echo -e "${Cwhite}$PROGNAME ${ARGS[@]} ${CNC}"
$PROGNAME "${ARGS[@]}"

res=$?
if [ $res -ne 0 ]
then
    echo "Error: can't execute command"
    exit 1;
fi

exit 0;
#}
