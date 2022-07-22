#!/bin/bash

USAGE="usage: $0 [-h] [--env ENV] [--] SCRIPT [ARGS]

activates the appropriate conda environment and calls python on SCRIPT
with ARGS

optional arguments:
  -h, --help  show this help message and exit
  --env ENV   Conda environment to activate before invoking SCRIPT
              (default: fipy)
  --log CONFIG LOG  Path to log configuration file template and
                    path for log file.
"

ENV=fipy

while [[ $# > 0 ]] && [[ $1 == -* ]]
do
    case "$1" in
        --env)
            ENV="$2"
            shift # option has parameter
            ;;
        --log)
            LOGCONFIG="$2"
            LOGFILE="$3"
            shift # option has two parameters
            shift
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

if [[ -n $LOGCONFIG ]]; then
    tmp_dir=$(mktemp -d -t fipylogconfig-XXXXXXXXXX)
    cp $LOGCONFIG $tmp_dir
    LOGCONFIG="${tmp_dir}/${LOGCONFIG##*/}"

    logpath=${LOGFILE%/*}
    mkdir -p $logpath
    if [[ -n $SLURM_JOB_ID ]]; then
        logbase=${LOGFILE##*/}
        logpref=${logbase%.*}
        logfext=${logbase##*.}

        LOGFILE="${logpath}/${logpref}.${SLURM_JOB_ID}.${logfext}"
    fi

    sed -i -e "s:%LOGFILE%:${LOGFILE}:g" $LOGCONFIG
fi

# https://stackoverflow.com/a/56155771/2019542
eval "$(conda shell.bash hook)"
conda activate $ENV
env FIPY_LOG_CONFIG=${LOGCONFIG} $@

if [[ -n $tmp_dir ]]; then
    rm -rf ${LOGCONFIG}
fi
