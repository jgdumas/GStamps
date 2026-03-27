#!/bin/bash
#########################################################################
# GStamps: C++ routines for the Local & the Global Postage Stamp Problem
# Authors:
#   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
#########################################################################

BEN=`dirname $0`
BIN=`echo ${BEN} | sed 's/benchmarks/bin/'`

PRGS="${BIN}/fibo ${BIN}/basis"

MAX=10
MIN=1
RAN=0

if [ "$#" -ge 1 ]; then
    MAX=$1
fi
if [ "$#" -ge 2 ]; then
    MIN=$2
fi
if [ "$#" -ge 3 ]; then
    RAN=$3
fi


# echo -n "k"
# for prg in ${PRGS}
#   do
#   echo -e -n '\t'`basename $prg`
# done
# echo

for j in $(seq 1 $MAX)
  do
  echo -n "$j"
  for prg in ${PRGS}
    do
    if [ "$RAN" -eq 0 ]; then
	echo -e -n "\t"`$prg $j $j 1|& grep nmax |cut -d' ' -f3`
    else
	echo -e -n "\t"`($prg $j $j| ${BIN}/range $j) 2> /dev/null`
    fi
  done
  echo
done
