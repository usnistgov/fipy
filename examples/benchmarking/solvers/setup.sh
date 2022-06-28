#!/bin/bash

USAGE="usage: $0 [-h] [--env ENV] [--] SCRIPT [ARGS]

activates the appropriate conda environment and calls python on SCRIPT
with ARGS

optional arguments:
  -h, --help  show this help message and exit
  --env ENV   Conda environment to activate before invoking SCRIPT
              (default: fipy)"

ENV=fipy

while [[ $# > 0 ]] && [[ $1 == -* ]]
do
    case "$1" in
        --env)
            ENV="$2"
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

if [[ "$#" < 1 ]]; then
    echo "$USAGE"
    exit 1
fi

source activate $ENV
$@
