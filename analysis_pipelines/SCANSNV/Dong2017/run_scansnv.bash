nice scansnv \
 --ref /net/refdata/HomoSapiens/hg19_CanonicalChr/genome_bwa-0.7.9a/hg19.fa \
 --dbsnp /net/refdata/HomoSapiens/hg19_GATKBundle-2.8/dbsnp_138.hg19.vcf \
 --shapeit-panel /net/refdata/HomoSapiens/hg19_shapeit/1000GP_Phase3 \
 --output-dir variant_calling/hg19/Clones_single-cells/sensitivity_1/ \
 --bam Clones ../Dong2017/mapping/hg19/samples/Clones.bps.sorted.bam \
 --bam IL-11 ../Dong2017/mapping/hg19/samples/IL-11.bps.sorted.bam \
 --bam IL-12 ../Dong2017/mapping/hg19/samples/IL-12.bps.sorted.bam \
 --bulk-sample Clones \
 --sc-sample IL-11 \
 --sc-sample IL-12 \
 --target-fdr 1.0 \
 --min-sc-alt 1 \
 --min-sc-dp 1 \
 --min-bulk-dp 1 \
 --regions-file gatk_regions_example.txt\
 --overwrite \
 --joblimit 44 --resume
