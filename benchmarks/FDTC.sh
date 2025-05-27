#!/bin/bash
#########################################################################
# GStamps: C++ routines for the Global Postage Stamp Problem
# Authors:
#   L. Colisson, J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
#########################################################################


# COVE=./bin/scover
# REAC=./bin/cover
REAC=./bin/reach
COVE=./bin/kcover
MOSS=./bin/scover
BASI=./bin/basis


MAXKpS=15
if [ "$#" -ge 1 ]; then
    MAXKpS=$1
fi

KBAS=35
if [ "$#" -ge 2 ]; then
    KBAS=$2
fi

TDoS=1
if [ "$#" -ge 3 ]; then
    TDoS=$3
fi


NUMI=7
if [ "$#" -ge 4 ]; then
    NUMI=$4
fi

JSTR=1
if [ "$#" -ge 5 ]; then
    JSTR=$5
fi


export LC_NUMERIC="en_US.UTF-8" # decimal with '.' 

function Median() {
    sorted_data=($(printf '%.6f\n' "${@}" | sort -n))
    local num_elements=${#sorted_data[@]}
    if (( $num_elements % 2 == 1 )); then     # Odd number of elements
        ((middle=$num_elements/2))
        val="${sorted_data[ $(($num_elements/2)) ]}"
    else                            # Even number of elements
        ((before_middle=$num_elements/2 - 1))
        ((after_middle=$num_elements/2))
        val=$(echo "scale=6;(${sorted_data[$before_middle]} + ${sorted_data[$after_middle]})/2" | bc -l)
    fi
    echo $val
}

function Run () {
    local vals=()
    for i in $(seq 1 $NUMI)
      do
      vals+=(`($1 $2 $3 0 | $4 $5 1 ) |& grep Cover | cut -d' ' -f3|sed 's/s//'`)
    done
    local res=$(Median ${vals[*]})
    echo -n ' '$res
}
  

if [ "$TDoS" -gt 0 ]; then 
    echo "Denominations=$KBAS, median $NUMI, from $JSTR to $MAXKpS" 1>&2 
    echo "s $REAC $COVE $MOSS" 1>&2 

    for j in $(seq $JSTR $MAXKpS)
      do
      JBAS=10
      if [ "$JBAS" -gt "$j" ]; then
	  JBAS=$j
      fi
      echo -n $j
      Run $BASI $KBAS $JBAS $REAC $j 
      Run $BASI $KBAS $JBAS $COVE $j 
      Run $BASI $KBAS $JBAS $MOSS $j 
      echo
    done

else
    echo "Stamps=$MAXKpS, median $NUMI, from $JSTR to $KBAS" 1>&2 
    echo "k $REAC $COVE $MOSS" 1>&2 

    for j in $(seq $JSTR $KBAS)
      do
      echo -n $j
      Run $BASI $j $MAXKpS $REAC $MAXKpS 
      Run $BASI $j $MAXKpS $COVE $MAXKpS 
      Run $BASI $j $MAXKpS $MOSS $MAXKpS 
      echo
    done
   
fi
