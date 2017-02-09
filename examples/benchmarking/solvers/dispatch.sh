#!/bin/bash

USAGE="usage: $0 [-h] [--env ENV] [--cmd CMD] [--np NP] [--] SCRIPT

Iterates over solvers and mesh sizes by calling setup.sh, which activates 
the appropriate conda environment and calls python on SCRIPT

positional arguments:
  SCRIPT      Python script to launch (expected to be in same
              directory as $0)

optional arguments:
  --cmd CMD   Command used to invoke SCRIPT, e.g., 'qsub -cwd' for
              Sun grid engine (default: bash)
  --env ENV   Conda environment to activate before invoking SCRIPT
              (default: fipy)
  --np NP     Number of processes to invoke SCRIPT with (default: 1)
  -h, --help  show this help message and exit"

CMD=bash
ENV=fipy
NP=1

while [[ $# > 0 ]] && [[ $1 == -* ]]
do
    case "$1" in
        --cmd)
            CMD="$2"
            shift # option has parameter
            ;;
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

if [[ "$#" != 1 ]]; then
    echo "$USAGE"
    exit 1
fi

SCRIPT=$1

for solver in trilinos scipy pysparse
do
    for size in 100 1000 10000 100000 1000000 # 10000000
    do
	${CMD} "${BASH_SOURCE%/*}/setup.sh" --env "${ENV}" --np "${NP}" -- $SCRIPT "--${solver}" "--numberOfElements=${size}"
    done
done
