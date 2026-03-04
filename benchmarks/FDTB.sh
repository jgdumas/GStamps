#!/bin/bash
#########################################################################
# GStamps: C++ routines for the Local & the Global Postage Stamp Problem
# Authors:
#   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
#########################################################################


PRGS="./bin/fibo ./bin/geom ./bin/alba ./bin/greedy ./bin/basis ./bin/dynprg"
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
  echo -e -n '\t'`basename $prg`
done
echo

for j in $(seq 1 $MAXKpS)
  do
  for k in $(seq $MINK $(( $j - 1 )))
    do
    s=$(( $j - $k ))
    if [[ $s -ge $MINS ]]
	then
	echo -n "$k $s"
	for prg in ${PRGS}
	  do
	  echo -e -n "\t"`$prg $k $s 1|& grep nmax |cut -d' ' -f3`
	done
	echo
    fi
  done
done
