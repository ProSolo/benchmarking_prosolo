#!/usr/bin/env bash

# provide all necessary files on the command line using globbing, e.g.:
# *.bwamem_UniVec.sorted.idxstats.txt

COLLECTION=collection.${1##*bwamem_}
OUTFILE=summary.${1##*bwamem_}

if [ -e $COLLECTION ]
then
  rm $COLLECTION
fi

for u in $@; do echo $u >>$COLLECTION; cat $u | awk '{ if ($3 > 10000) { print $0 } }' >>$COLLECTION; done

echo -e "#VectorName\tVectorLength\tTotVectorCov\tAvgVectorCovInCovSmps\tNumCovSmps\tCovSmps" >$OUTFILE
cat $COLLECTION | awk ' 
  BEGIN { OFS="\t" } 
  { 
    if ($1 ~ /sorted\.idxstats\.txt/) 
      { SMP=$1 } 
    else
      { V_COUNT[$1]++; V_SMP[$1] = V_SMP[$1] "," SMP; V_LEN[$1] = $2; V_COV[$1] += $3 }
  }
  END {
    for (v in V_COUNT)
      { print v,V_LEN[v],V_COV[v],V_COV[v]/V_COUNT[v],V_COUNT[v],V_SMP[v] }
  }' >>$OUTFILE
