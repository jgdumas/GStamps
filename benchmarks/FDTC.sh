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


MAXKpS=10	# K or Kmax: Denominations
if [ "$#" -ge 1 ]; then
    MAXKpS=$1
fi

KBAS=15		# S or Smax: Stamps
if [ "$#" -ge 2 ]; then
    KBAS=$2
fi

TDoS=1		# Loop over Denominations or Stamps
if [ "$#" -ge 3 ]; then
    TDoS=$3
fi


NUMI=7		# Number of iterations from which take the median
if [ "$#" -ge 4 ]; then
    NUMI=$4
fi

JSTR=1		# Start of loop (Kmin or Smin)
if [ "$#" -ge 5 ]; then
    JSTR=$5
fi

PREC=6		# Precision of the reported timings
if [ "$#" -ge 6 ]; then
    PREC=$6
fi

export LC_NUMERIC="en_US.UTF-8" # decimal with '.'

function Median() {
    local sorted_data=($(printf '%f\n' "${@}" | sort -n))
    local num_elements=${#sorted_data[@]}
    if (( $num_elements % 2 == 1 )); then     # Odd number of elements
	((middle=$num_elements/2))
	local val="${sorted_data[ $(($num_elements/2)) ]}"
    else                            # Even number of elements
	((before_middle=$num_elements/2 - 1))
	((after_middle=$num_elements/2))
	local val=$(echo "scale=${PREC};(${sorted_data[$before_middle]} + ${sorted_data[$after_middle]})/2" | bc -l)
    fi
    printf '%.*f' $PREC $val
}


function count_unique() {
  local LC_ALL=C

  if [ "$#" -eq 0 ] ; then
    echo 0
  else
    echo "$(printf "%s\000" "$@" |
	      sort --zero-terminated --unique |
	      grep --null-data --count .)"
  fi
}

GRE='\033[0;32m'
RED='\033[0;41m'
NC='\033[0m'    # No Color

function Check () {
    local vals=()
    for prg in $5 $6 $7
      do
      vals+=(`($1 $2 $3 0 | $prg $4 1) 2> /dev/null`)
    done
    if [ "$(count_unique "${vals[@]}")" -eq 1 ] ; then
	echo -e -n "${GRE}$4 $2 success${NC} "
    else
	echo -e -n "${RED}$4 $2 ERRORS:***${NC} ${vals[@]} ${RED}***${NC} "
    fi
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
    echo -e -n "${GRE}k s correctness${NC} " 1>&2
    echo "s $REAC $COVE $MOSS" 1>&2
    (printf '=%.0s' {1..50} ; printf '\n') 1>&2


    for j in $(seq $JSTR $MAXKpS)
      do
      JBAS=10
      if [ "$JBAS" -gt "$j" ]; then
	  JBAS=$j
      fi
      Check $BASI $KBAS $JBAS $j $REAC $COVE $MOSS 1>&2
      echo -n $j
      Run $BASI $KBAS $JBAS $REAC $j
      Run $BASI $KBAS $JBAS $COVE $j
      Run $BASI $KBAS $JBAS $MOSS $j
      echo
    done

else
    echo "Stamps=$MAXKpS, median $NUMI, from $JSTR to $KBAS" 1>&2
    echo -e -n "${GRE}k s correctness${NC} " 1>&2
    echo "k $REAC $COVE $MOSS" 1>&2
    (printf '=%.0s' {1..50} ; printf '\n') 1>&2

    for j in $(seq $JSTR $KBAS)
      do
      Check $BASI $j $MAXKpS $MAXKpS $REAC $COVE $MOSS 1>&2
      echo -n $j
      Run $BASI $j $MAXKpS $REAC $MAXKpS
      Run $BASI $j $MAXKpS $COVE $MAXKpS
      Run $BASI $j $MAXKpS $MOSS $MAXKpS
      echo
    done

fi

echo 1>&2
