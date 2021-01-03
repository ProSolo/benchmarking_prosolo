#!/usr/bin/env bash

SNAKEMAKE_DIR="../"
# directory where accessory scripts for pipeline rules live
SCRIPTS_DIR="$SNAKEMAKE_DIR/scripts"
# make scripts available
PATH="$SCRIPTS_DIR:$PATH"
# create a conda environment for your pipeline with all the necessary tools installed
eval "$( conda shell.posix activate Benchmarking_ProSolo )"

RUNS=( `ls *stderr.log` )
STDERR=${#RUNS[@]}stderr.log
STATS=${#RUNS[@]}stats.json
ENV=${#RUNS[@]}conda-package-list.txt

conda list --export >$ENV

#export GRB_LICENSE_FILE=/path/to/gurobi_license/gurobi.lic
LOGS=sge_logs/${#RUNS[@]}/

#mkdir -p $LOGS

##snakemake --cluster "qsub -V -cwd -pe multislot {threads} -q all.q -o ./$LOGS -e ./$LOGS" --use-conda --output-wait 61 --timestamp --stats $STATS --cores 180 2> >(tee $STDERR >&2)
#
### alternatively use the below command if you have rules that specify a per-thread memory usage in GBs in resources.mem_gb
#snakemake --cluster "qsub -V -cwd -pe smp {threads} -l h_vmem={resources.mem_gb}G -q compute.q -o ./$LOGS -e ./$LOGS" --use-conda --conda-frontend mamba --output-wait 68 --stats $STATS --cores 60 2> >(tee $STDERR >&2)
nice snakemake --use-conda --conda-frontend mamba --resources mem_gb=300 --output-wait 31 --stats $STATS --cores 40 2> >(tee $STDERR >&2)

## alternatively use the below command for detailed debugging
#snakemake --debug --timestamp 2> >(tee $STDERR >&2)
eval "$( conda shell.posix deactivate )"
