#!/usr/bin/env bash

SNAKEMAKE_DIR="../"
## load all needed tools
# directory where accessory scripts for pipeline rules live
SCRIPTS_DIR="$SNAKEMAKE_DIR/scripts"
# source snakemake .tools-files to make executables that you cannot install via conda available as required by snakemake rules
source script/scripts.tools
# create a conda environment for your pipeline with all the necessary tools installed
source activate Benchmarking_ProSolo

RUNS=( `ls *stderr.log` )
STDERR=${#RUNS[@]}stderr.log
STATS=${#RUNS[@]}stats.json
ENV=${#RUNS[@]}conda-package-list.txt

conda list --export >$ENV

LOGS=sge_logs/${#RUNS[@]}/

mkdir -p $LOGS

snakemake --cluster "qsub -V -cwd -pe multislot {threads} -q all.q -o ./$LOGS -e ./$LOGS" --use-conda --output-wait 41 --timestamp --stats $STATS --cores 160 2> >(tee $STDERR >&2)

## alternatively use the below command for detailed debugging
#snakemake --debug --timestamp 2> >(tee $STDERR >&2)
source deactivate
