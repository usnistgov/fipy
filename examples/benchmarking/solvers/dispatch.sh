#!/bin/bash

for solver in trilinos scipy pysparse
do
    for size in 100 1000 10000 100000 1000000 10000000
    do
	qsub -cwd script.sh "--$solver" --numberOfElements=$size
    done
done
