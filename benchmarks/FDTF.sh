#!/bin/bash
#########################################################################
# GStamps: C++ routines for the Local & the Global Postage Stamp Problem
# Authors:
#   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
#########################################################################

BEN=`dirname $0`
BIN=`echo ${BEN} | sed 's/benchmarks/bin/'`

PRGS="${BIN}/basis"
DPTH="${BIN}/depthrange"
DECO="${BIN}/decompositions"

MAXK=5
STPS=4

if [ "$#" -ge 1 ]; then
    MAXK=$1
fi
if [ "$#" -ge 2 ]; then
    STPS=$2
fi

function is_power_of_two() {
    declare -i n=$1
    (( n > 0 && (n & (n - 1)) == 0 ))
}
function is_power_of_two_mOne() {
    declare -i n=$1
    (( (n+1)> 0 && (n & (n + 1)) == 0 ))
}

GRE='\033[0;32m'
RED='\033[0;41m'
NC='\033[0m'    # No Color

function test_k_s() {
    _k=$1
    _s=$2

    depth=(`(./bin/basis ${_k} ${_s} | ./bin/depthrange ${_s}) 2> /dev/null`)
#     echo " ----------- "
#     echo "${#depth[@]}: ${depth[@]}"

    decom=(`(./bin/basis ${_k} ${_s} | ./bin/decompositions ${_s} | sed 's/[^=+]//g' | sort | uniq -c| sed -r 's/[^0-9]+//g' | tr '\n' ' ') 2> /dev/null`)
#      echo "${#decom[@]}: ${decom[@]}"
#      echo " ----------- "

    sum=0
    ide=0
    errors=0
    for i in ${!decom[@]}; do
	((sum += ${decom[i]}))
	if is_power_of_two_mOne $i; then
	    if [[ $sum -eq ${depth[ide]} ]]; then
		echo -en "${GRE}$sum${NC} "
	    else
		((errors++))
		echo -en "${RED}$sum != ${depth[ide]}${NC} "
	    fi
	    sum=0
	    ((ide++))
	fi
    done
    if ! is_power_of_two ${#decom[@]}; then
	    if [[ $sum -eq ${depth[-1]} ]]; then
		echo -en "${GRE}$sum${NC} "
	    else
		((errors++))
		echo -en "${RED}$sum != ${depth[-1]}${NC} "
	    fi
    fi
    if [[ $errors -gt 0 ]]; then
	echo -e "${RED}*** ERROR ***${NC}"
    else
	echo -e "${GRE}SUCCESS${NC}"
    fi
}

for k in $(seq 1 $MAXK)
  do
  for s in $(seq 1 $STPS)
    do
    echo -n "<$k,$s>: "
    test_k_s $k $s
  done
  echo
done
