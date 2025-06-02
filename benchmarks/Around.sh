#!/bin/bash
#########################################################################
# GStamps: C++ routines for the Global Postage Stamp Problem
# Authors:
#   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
#########################################################################

K=$1
S=$2

G=2
if [ "$#" -ge 3 ]; then
    G=$3
fi

SUPP=./bin/supplement
COMP=./bin/complement
DYNP=./bin/basis
RESF=/tmp/res_Around.$$

##### 
echo "## $0 $*: $RESF"


if [ "$#" -ge 4 ]; then
    read -r BASIS
else
    BASIS=`$DYNP $K $S 1`
fi

echo "## Base: $BASIS"

for j in $(seq -$G $G)
  do
  for i in $(seq -$G $G)
    do
    ik=$(( $K + $i ))
    js=$(( $S + $j ))
    if [ "${ik}" != ${K} ] || [ "${js}" != "${S}" ]; then
#####
	mk=$(( $ik - $G ))
	echo "#======================================#"| tee -a $RESF
	CMD="$COMP $ik $js 1"
	SEQ=`echo ${BASIS} | awk -v max=${mk} '{out=$1; for(i=2;i<=max;i++) {out=out" "$i}; print out}'`
	echo "# $ik $js $mk : $SEQ | $CMD" |& tee -a $RESF
	FOUND=`(echo "$SEQ" | $CMD ) 2>> $RESF`
	echo "## Found: `grep nmax $RESF|tail -1|cut -d' ' -f4` (${FOUND})"
    fi
  done
done

