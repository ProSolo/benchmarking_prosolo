#!/usr/bin/env bash

# provide all the necessary files on the command line using globbing, e.g.:
# *.bwamem_UniVec.stats.txt

OUTFILE=summary.${1##*bwamem_}

if [ -e $OUTFILE ]
then
  rm $OUTFILE
fi

for u in $@; do echo $u >>$OUTFILE; cat $u | grep ^SN | cut -f 2- | grep -P "^reads (un)?mapped:" >>$OUTFILE; done
