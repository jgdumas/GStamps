#!/bin/bash
#########################################################################
# GStamps: C++ routines for the Local & the Global Postage Stamp Problem
# Authors:
#   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
#########################################################################
# Example:
# > ./benchmarks/format.sh ./bin/basis 16 6
# Example (with /tmp/t echoing 1 6 7):
# > ./benchmarks/format.sh /tmp/t 3 5


RES=(`$* | cut -d' ' -f1-$2`)
MAX=`echo ${RES[@]} | ./bin/cover $3`

echo "${#RES[@]} $3 $MAX "`echo "${RES[@]}"|sed 's/ /,/g'` | ./benchmarks/awk.gb
