#!/bin/env bash

print_help() {
    print_usage
    cat <<HELP_TEXT
Initialize the set of scripts by setting appropriate 
PATH and PYTHONPATH variablies.

positional arguments:
  NAME                  name of the set

options:
  -h, --help            show this help message and exit
  -u, --usage           show usage information
  -l, --list            show available sets

Also, you can pass 'all' to initialize all the available 
scripts or 'python' to set only the PYTHONPATH. If no 
arguments are passed, the 'all' option will be assumed.
HELP_TEXT
}

print_usage() {
echo Usage:
echo "   source init.sh [OPTIONS] NAME"
}

print_targers() {
    echo "Available sets:"
    echo "  ${targets[@]}"
}

append_path() {
    # $1 - name of variable to modify (e.g. PATH)
    # $2 - string to append
    orig=${!1}
    if [[ "$orig" =~ $2 ]] #the string is already present, nothing to do
    then
        echo $orig  
    else    #the string is absent, append
        echo "$orig:$2"
    fi
}

init_collections() {
    allowed_keys=(${targets[*]} "python")
    for key in $@   #Test the key is in allowed
    do
        if [[ ! " ${allowed_keys[*]} " =~ " $key " ]]
        then
            echo "Error: the set named '${key}' is not found."
            return 1 2>/dev/null || exit 1
        fi
    done
    for key in $@  #Apply
    do
        [ $key == "python" ] && continue
        export PATH=$(append_path "PATH" "${SCRIPT_DIR}/libexec/$key")
    done
}

#### main ####

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)
targets=($(ls -1 $SCRIPT_DIR/libexec))

[ $# -lt 1 ] && set -- all

for opt in $(getopt -u -o "h,u,l" --longoptions "help,usage,list" -n $0 -- "$@")
do
    case $opt in
        all)
            set -- ${targets[@]}
            break
            ;;
        -h | --help)
            print_help
            return 0 2>/dev/null || exit 0
            ;;
        -l | --list)
            print_targers
            return 0 2>/dev/null || exit 0
            ;;
        -u |--usage)
            print_usage
            return 0 2>/dev/null || exit 0
    esac
done

init_collections $@
export PYTHONPATH=$(append_path "PYTHONPATH" "${SCRIPT_DIR}/python-packages")

