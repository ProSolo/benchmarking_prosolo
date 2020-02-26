#!/usr/bin/env bash

# make conda work
source /net/cancergenomics/miniconda3/etc/profile.d/conda.sh

SNAKEMAKE_DIR="../"
# make scripts available
PATH="$SCRIPTS_DIR:$PATH"
# create a conda environment for your pipeline with all the necessary tools installed
conda activate Benchmarking_ProSolo

RUNS=( `ls *stderr.log` )
STDERR=${#RUNS[@]}stderr.log
STATS=${#RUNS[@]}stats.json
ENV=${#RUNS[@]}conda-package-list.txt

conda list --export >$ENV

#LOGS=sge_logs/${#RUNS[@]}/
#
#mkdir -p $LOGS

##snakemake --cluster "qsub -V -cwd -pe multislot {threads} -q all.q -o ./$LOGS -e ./$LOGS" --use-conda --output-wait 61 --timestamp --stats $STATS --cores 180 2> >(tee $STDERR >&2)
#
### alternatively use the below command if you have rules that specify a per-thread memory usage in GBs in resources.mem_gb
#snakemake --cluster "qsub -V -cwd -pe smp {threads} -l h_vmem={resources.mem_gb}G -q compute.q -o ./$LOGS -e ./$LOGS" --output-wait 31 --stats $STATS --cores 108 2> >(tee $STDERR >&2)
nice snakemake --use-conda --output-wait 31 --stats $STATS --cores 40 2> >(tee $STDERR >&2)

## alternatively use the below command for detailed debugging
#snakemake --debug --timestamp 2> >(tee $STDERR >&2)
conda deactivate
