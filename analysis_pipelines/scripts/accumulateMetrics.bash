#!/usr/bin/env bash

# arguments for this script:
# $1 - file suffix for files to accumulate
# $2 - lines starting with this Perl regex will be kept only from the first file,
#      i.e. these should be the starting characters of all header lines except
#      those starting with a '#'
# $3 and following - list of input files

SUFF=$1
shift
HEAD=$1
shift

FILES=( $@ )

cat ${FILES[0]} | grep -P "^(#|${HEAD})" 

for f in ${FILES[@]}
do
    run=${f##*/}
    run=${run%%$SUFF}
    cat $f | grep -P -v "^(#|$|\t0\t0|${HEAD})" | sed -e "s/^[^\t]*/$run/" 
done

