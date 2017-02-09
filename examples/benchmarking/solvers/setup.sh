#!/bin/bash

USAGE="usage: $0 [-h] [--env ENV] [--np NP] [--] SCRIPT [ARGS]

activates the appropriate conda environment and calls python on SCRIPT

positional arguments:
  SCRIPT      Python script to launch (expected to be in same
              directory as $0)
  ARGS        Arguments to pass to SCRIPT

optional arguments:
  --env ENV   Conda environment to activate before invoking SCRIPT
              (default: fipy)
  --np NP     Number of processes to invoke SCRIPT with (default: 1)
  -h, --help  show this help message and exit"

NP=1
ENV=fipy

while [[ $# > 0 ]] && [[ $1 == -* ]]
do
    case "$1" in
        --env)
            ENV="$2"
            shift # option has parameter
            ;;
        --np)
            NP="$2"
            shift # option has parameter
            ;;
        -h|--help)
            echo "$USAGE"
            exit 0
            ;;
        --)
            # end of options
            shift
            break
            ;;
        -*)
            # unknown option
            echo Unknown option: $1>&2
            exit 10 
            ;;
    esac
    shift # option(s) fully processed, proceed to next input argument
done

echo "$@"

if [[ "$#" < 1 ]]; then
    echo "$USAGE"
    exit 1
fi

SCRIPT=$1
shift

if [[ $NP > 1 ]]; then
    MPI="mpirun -np ${NP}"
else
    MPI=""
fi

source activate $ENV
${MPI} python "${BASH_SOURCE%/*}/${SCRIPT}" "$@"
