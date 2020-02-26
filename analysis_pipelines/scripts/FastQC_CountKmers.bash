#! /usr/bin/env bash

echo "## Kmer - k-mer sequence"
echo "## #Samples - number of samples the k-mer was found in"
echo "## Tot#Reads - total number of reads the k-mer was found in"
echo "## Avg#Reads - average number of reads the k-mer was found across the samples it was found in"
echo "## AvgOverrepAt - average ratio of over-representation at maximum over-representation read positions"
echo "## MaxPos - maximum over-representation read positions in all the samples it was found in"
echo "## "
echo "## Kmer #Samples Tot#Reads Avg#Reads AvgOverrepAt MaxPos"

xmllint --html $@ | grep -A 4 -P "<td>[ACGT]{7}</td>" | awk 'BEGIN { RS="\n--\n"; FS="(</td>\n<td>|<td>|</td>)" } { print $2,$3,$4,$5,$6 }' | awk '{ smp_count[$1]++; tot_count[$1]+=$2; max_obs_exp[$1]+=$4; if (smp_count[$1] == 1) { pos[$1] = $5 } else { pos[$1] = (pos[$1] "," $5) } } END { for (c in smp_count) { print c, smp_count[c], tot_count[c], tot_count[c]/smp_count[c], max_obs_exp[c]/smp_count[c], pos[c] } }' | sort -n -k2,2
