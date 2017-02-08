#!/bin/bash

source activate condafipy
python "${BASH_SOURCE%/*}/solvers.py" "$@"
