# Snakemake pipelines for the benchmarking of ProSolo

This folder contains the Snakemake pipelines used for the ProSolo benchmarking paper "ProSolo:
Accurate Variant Calling from Single Cell DNA Sequencing Data".

DISCLAIMER: These analysis pipelines have grown over the years, and I cannot recommend them
as a model for new pipelines. For some ready-to-use workflows and good templates for new
workflows, see the [Snakemake Workflow project](https://github.com/snakemake-workflows/docs).


## General setup

To run any of this, you will have to have a (bio-)conda installation present. I recommend
following the [bioconda installation instructions](https://bioconda.github.io/user/install.html)
and then adding in mamba for faster dependency resolution:

    conda install mamba

Once this is set up, you need to create an environment that contains `snakemake` and some other
(older) dependencies, using the environment specification provided in this folder:

    mamba env create --file Benchmarking_ProSolo_minimal.yml

(An older version of this specification, which pins all packages and thus makes dependency
resolution harder or even impossible, is also available for documentation:
`Benchmarking_ProSolo.explicit.txt`.)

Once you have this environment, the four datasets should be run through their respective
pipelines in the following order:

1. `Dong2017`: Two whole genome sequenced single cells and a corresponding bulk, with the ground
   truth generation within the same workflow.
2. `Hoell2014`: The pedigree dataset from which to generate the ground truth for `Laehnemann2017`.
3. `Laehnemann2017`: Five whole exome sequenced single cells and a corresponding bulk.
4. `Wang2014`: 16 normal single cells and 16 tumor single cells with two corresponding bulk
   sequencing samples, with ground truth clonal and subclonal variants provided through targeted
   resequencing validation in the original publication (Supplementary Tables 6 and 7).

For each dataset (`Dong2017`, `Hoell2014`, `Laehnemann2017`, `Wang2014`), a separate `Snakefile`,
`config.yaml` and `pipeline.wrapper.bash` exist in the respective subfolders. You will have to adjust:

### `config.yaml`

Generally, all the paths have to be set correctly for your respective environment. Search for
`/path/to/` and adjust accordingly. For example, a `tmp`-path should be set to a location with
lots of usable storage for some larger intermediate outputs.

Further, adjust paths for the following:

#### dataset

You will have to download the data for each pipeline, put it into a subdirectory of the
pipeline folder (e.g. `Dong2017/data/`) and set the `folders`->`data` path accordingly. Namely,
you will have to get:

1. `Dong2017`: From [SRA BioProject PRJNA305211](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA305211&ff=on)
   download samples `IL-1c`, `IL-11`, `IL-12`, `Clone-1`, `Clone-2` and `Clone-3` (they have
   the accessions `SRR2976561, SRR2976562, SRR2976563, SRR2976564, SRR2976565, SRR2976566`).
2. `Hoell2014`: From [EGA Dataset EGAD00001005929](https://www.ebi.ac.uk/ega/datasets/EGAD00001005929)
   download samples `EGAN00002446901, EGAN00002446902, EGAN00002446903, EGAN00002446904, EGAN00002446905, EGAN00002446906`.
3. `Laehnemann2017`: From the same [EGA Dataset EGAD00001005929](https://www.ebi.ac.uk/ega/datasets/EGAD00001005929)
   further download samples `EGAN00002446895, EGAN00002446896, EGAN00002446897, EGAN00002446898, EGAN00002446899, EGAN00002446900`.
4. `Wang2014`: From [SRA Bioproject PRJNA168068](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA168068&f=sample_name_s%3An%3Atnbc%3Blibraryselection_s%3An%3Amda%2Cpcr%3Ac&o=acc_s%3Aa&s=SRR1163035,SRR1163012,SRR1163013,SRR1163019,SRR1163026,SRR1163027,SRR1163034,SRR1163043,SRR1163053,SRR1163070,SRR1163074,SRR1163083,SRR1163084,SRR1163091,SRR1163095,SRR1163148,SRR1163149,SRR1163150,SRR1163151,SRR1163152,SRR1163153,SRR1163154,SRR1163155,SRR1163156,SRR1163157,SRR1163158,SRR1163159,SRR1163160,SRR1163161,SRR1163162,SRR1163163,SRR1163164,SRR1298936,SRR1163508#)
   download samples `SRR1163035,SRR1163012,SRR1163013,SRR1163019,SRR1163026,SRR1163027,SRR1163034,SRR1163043,SRR1163053,SRR1163070,SRR1163074,SRR1163083,SRR1163084,SRR1163091,SRR1163095,SRR1163148,SRR1163149,SRR1163150,SRR1163151,SRR1163152,SRR1163153,SRR1163154,SRR1163155,SRR1163156,SRR1163157,SRR1163158,SRR1163159,SRR1163160,SRR1163161,SRR1163162,SRR1163163,SRR1163164,SRR1298936,SRR1163508`.

#### references

You will have to provide all the respective reference files somewhere on your file system and
set the paths accordingly (search for `/path/to/refdata/`).


The biggest amount is provided by the
[GATK Resource Bundle](https://gatk.broadinstitute.org/hc/en-us/articles/360036212652-Resource-Bundle).
Follow the instructions in the linked documentation to download `b37/hg19` for comparable results.

For the `Hoell2014` and the `Laehnemann2017` pipelines, you will also need to download the
`Agilent SureSelect SureSelectXT Human All Exon V5+UTRs Kit` by
[creating a SureDesign account](https://earray.chem.agilent.com/suredesign/help/Set_up_an_account.htm)
and following the instructions to
[download design files](https://earray.chem.agilent.com/suredesign/help/WebHelp.htm#Target_enrichment_design_files_available_for_download.htm)
to get `S04380219_Covered.bed`. You will have to manually turn this into an interval list file using
a [picardtools](https://broadinstitute.github.io/picard/) generated dict file of the previously
downloaded reference genome:

    conda create -n picard picard
    conda activate picard
    picard CreateSequenceDictionary REFERENCE=/path/to/refdata/hg19.fa OUTPUT=REFERENCE=/path/to/refdata/hg19.dict
    picard BedToIntervalList INPUT=S04380219_Covered.bed SEQUENCE_DICTIONARY=/path/to/refdata/hg19.dict OUTPUT=SureSelectHumanAllExonV5plusUTR.target.hg19.interval_list
    grep -v "^@" SureSelectHumanAllExonV5plusUTR.target.hg19.interval_list > SureSelectHumanAllExonV5plusUTR.target.noHeader.hg19.interval_list

Point the `settings`->`capture` entry in the respective `config.yaml`s to the location containing
the two above created versions of the file.

The same has to be done with a different capture kit for the `Wang2014` dataset.
The respective BED file is no longer provided by Illumina, but can be found [via this Biostars answer](https://www.biostars.org/p/144554/#144561).

### `pipeline.wrapper.bash`

You will have to adjust the `snakemake` call to your local system. The versions provided here
were created for different server environments. `Hoell2014/pipeline.wrapper.bash` is an example
for using snakemake with an SGE grid engine using the `qsub` command. The other two pipelines
also contain the SGE commands (commented out), but have an active version of the command
that can be executed directly on a machine with the respective amount of `--cores` available.
When editing this command, make sure that the `--use-conda` argument remains part of it.

### Changing the pipeline

If you want to change the pipeline, or have to debug things that don't work with the above
described setup, all of what the pipeline does lives in the `rules` directory. It contains all
the Snakemake rules necessary for all three pipelines. 

## Reproducing SCAN-SNV results

SCAN-SNV is distributed via conda, but not as a stand-alone tool, but as a wrapper around a
Snakemake pipeline. See the
[SCAN-SNV documentation](https://github.com/parklab/scan-snv/blob/6df4c19d6cd3cd0b796bf858eea874f79df34ad3/README.md)
for details. This setup made it impossible to integrate SCAN-SNV into the existing
pipelines and SCAN-SNV had to be installed into its own environment and run separately. For
the installation, please follow the instructions linked above and possibly modify them by
considering the info provided in [SCAN-SNV issue #4](https://github.com/parklab/scan-snv/issues/4).

For the exact `config.yaml` and `Snakefile` used, see the subdirectories of `SCAN-SNV/`. As
we tried to increase SCAN-SNV's sensitivity, we had to change the `Snakefile`. For `scansnv`
to work as we ran it, this requires that you move one of the (identical) `Snakefile` (from
`SCANSNV/Dong2017` or `SCANSNV/Laehnemann2017`) into the `scansnv` environment created while
following the above install instructions. In particular, you will have to:

    cp SCANSNV/Dong2017/Snakefile /path/to/miniconda3/envs/scansnv/lib/scansnv/

Sorry for the inconvenience -- this is the only way we found to make this work. Once this is
in place, check `SCAN-SNV/Dong2017/run_scansnv.bash` for the exact command that was run for
that dataset and run the `Laehneman2017` datasets accordingly. Please note, that
the respective main pipeline for that dataset will have to have run up to the point where the
necessary BAM files have been generated. Please also note, that the recent installation
problems seem to have created dependency issues for using `r-tibble` in the analysis of
SCAN-SNV results. We were thus unable to process SCAN-SNV results systematically for the
`Wang2014` dataset and excluded it from this analysis. 
