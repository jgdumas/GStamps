#!/bin/bash
#########################################################################
# GStamps: C++ routines for the Local & the Global Postage Stamp Problem
# Authors:
#   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
#########################################################################

BEN=`dirname $0`
BIN=`echo ${BEN} | sed 's/benchmarks/bin/'`

PRGS="${BIN}/bala ${BIN}/basis ${BIN}/brute"
MAXKpS=10

MINS=2
MINK=2


if [ "$#" -ge 1 ]; then
    MAXKpS=$1
fi
if [ "$#" -ge 2 ]; then
    MINK=$2
fi
if [ "$#" -ge 3 ]; then
    MINS=$3
fi


echo -n "k s"
for prg in ${PRGS}
  do
  echo -n " `basename $prg`"
done
echo

GRE='\033[0;32m'
RED='\033[0;41m'
NC='\033[0m'    # No Color

for j in $(seq 1 $MAXKpS)
  do
  for k in $(seq $MINK $(( $j - 1 )))
    do
    s=$(( $j - $k ))
    if [[ $s -ge $MINS ]]
	then
	echo -n "$k $s "
	rg=()
	for prg in ${PRGS}
	  do
	  rg+=(`$prg $k $s 1|& grep nmax |cut -d' ' -f3`)
	done
	echo -n ${rg[*]}
	if [[ ${rg[2]} -eq ${rg[1]} && ${rg[1]} -ge ${rg[0]} ]]
	    then
	    echo -e " ${GRE}success${NC}"
	    else
	    echo -e " ${RED} *** FAIL *** {NC}"
	fi
    fi
  done
done
