#!/bin/bash

RES=(`$* | cut -d' ' -f1-$2`)
MAX=`echo ${RES[@]} | ./cover $3`

echo "${#RES[@]} $3 $MAX "`echo "${RES[@]}"|sed 's/ /,/g'` | awk.gb
