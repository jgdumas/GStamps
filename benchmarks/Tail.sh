#!/bin/bash
####################################################################
# GStamps: C++ routines for the Global Postage Stamp Problem
# Authors: J-G. Dumas, A. Galan, B. Grenet, A. Maignan, D. S. Roche
####################################################################

K=$1
S=$2

G=2
if [ "$#" -ge 3 ]; then
    G=$3
fi

SUPP=./bin/supplement
COMP=./bin/complement
DYNP=./bin/dynprg
RESF=/tmp/res_Tail.$$

##### 
echo "## $0 $*: $RESF"


##### 
echo "#======================================#"| tee -a $RESF
CMD="$COMP $K $S 1 $G"
echo "## $0 $1 $2: `date`: (c ): $CMD" |& tee -a $RESF
$CMD &>> $RESF

##### 
echo "#======================================#"| tee -a $RESF
MAXS=$(( $S + 1 ))
for j in $(seq 2 $MAXS)
  do
  CMD="$SUPP $K $S 2 $G"
  echo "## $0 $1 $2: `date`: 1 $j: $CMD" |& tee -a $RESF
  (echo "1 $j" | $CMD) &>> $RESF

done


##### 
CMD="$COMP $K $S 1"
echo "#======================================#"| tee -a $RESF
for j in $(seq -$G $G)
  do
  jk=$(( $K + $j ))
  mk=$(( $K - $G ))
  SEQ=`$DYNP $jk $S 1 | awk -v max=$mk '{out=$1; for(i=2;i<=max;i++) {out=out" "$i}; print out}'`
  echo "# $jk $S $mk : $SEQ | $CMD" |& tee -a $RESF
  echo "## $0 $1 $2: `date`: k $j: $CMD 0" |& tee -a $RESF
  (echo "$SEQ" | $CMD ) &>> $RESF
done

##### 
echo "#======================================#"| tee -a $RESF
for j in $(seq 1 $G)
  do
  js=$(( $S - $j ))
  mk=$(( $K - $G ))
  SEQ=`$DYNP $K $js 1 | awk -v max=$mk '{out=$1; for(i=2;i<=max;i++) {out=out" "$i}; print out}'`
  echo "# $K $js $mk : $SEQ | $CMD" |& tee -a $RESF
  echo "## $0 $1 $2: `date`: $j s: $CMD 0" |& tee -a $RESF
  (echo "$SEQ" | $CMD ) &>> $RESF
done


echo "#======================================#"| tee -a $RESF
egrep '(##|nmax|^1 )' $RESF | sed 's/Cpmt /Cpmt/g' | paste -d " " - - - | awk '{out=$1" "$3" "$4" "$20" <-- "; for (i = 22; i <= NF; i++) {out = out" "$i}; print out}' | sort -k4n