#!/bin/bash
#########################################################################
# GStamps: C++ routines for the Global Postage Stamp Problem
# Authors:
#   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
#########################################################################

RES=(`$* | cut -d' ' -f1-$2`)
MAX=`echo ${RES[@]} | ./bin/cover $3`

echo "${#RES[@]} $3 $MAX "`echo "${RES[@]}"|sed 's/ /,/g'` | ./benchmarks/awk.gb
