#!/bin/bash

USAGE="usage: $0 [-h] [--env ENV] [--cmd CMD] [--np NP] [--mprof] [--] SCRIPT [ARGS]

Iterates over solvers and mesh sizes by calling setup.sh, which activates
the appropriate conda environment and calls python on SCRIPT

positional arguments:
  SCRIPT      Python script to launch (expected to be in same
              directory as $0)

optional arguments:
  -h, --help  show this help message and exit
  --qsub      Invoke SCRIPT using 'qsub -cwd' for Sun grid engine
              (default: invoke using bash)
  --env ENV   Conda environment to activate before invoking SCRIPT
              (default: fipy)
  --np NP     Number of processes to invoke SCRIPT with (default: 1)
  --mprof     Whether to run mprof profiler (default: False)"

QSUB=0
ENV=fipy
NP=1
PYTHON=python

while [[ $# > 0 ]] && [[ $1 == -* ]]
do
    case "$1" in
        --qsub)
            QSUB=1
            ;;
        --env)
            ENV="$2"
            shift # option has parameter
            ;;
        --np)
            NP="$2"
            shift # option has parameter
            ;;
        --mprof)
            PYTHON="mprof run"
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

for library in trilinos scipy pysparse
do
    for solver in cg pcg cgs gmres lu
    do
        for size in 100 1000 10000 100000 1000000 10000000
        do
            dir="Data/`uuidgen`"
            mkdir -p $dir
            INVOCATION="${MPI} ${PYTHON} ${BASH_SOURCE%/*}/${SCRIPT} \
              --${library} --numberOfElements=${size} --solver=${solver} --output $dir $@"
            if [[ $QSUB == 1 ]]; then
                qsub -cwd -pe nodal ${NP} -q "wide64" -o "${dir}" -e "${dir}" "${BASH_SOURCE%/*}/setup.sh" --env "${ENV}" -- ${INVOCATION}
            else
                bash "${BASH_SOURCE%/*}/setup.sh" --env "${ENV}" -- ${INVOCATION}
            fi
        done
    done
done
