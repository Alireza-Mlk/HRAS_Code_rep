# HRAS rMATS Analysis

Pipeline and [R](https://www.r-project.org/) codes used for HRAS splicing analysis.  
The pipeline is written using [CGAR-CORE](https://github.com/cgat-developers/cgat-core). It uses [rMATS](https://github.com/Xinglab/rmats-turbo) and analyses each given time point against all other 
time points in a series of two-group analyses.  
Requires [STAR](https://github.com/alexdobin/STAR) aligned BAM files with the naming scheme, [SAMPLE NAME]-T[SAMPLE TIME]R-[SAMPLE REPLICATE].star.bam (e.g, HRAS-T1-R1.star.bam), and a GTF file to run.  
The GSEA analysis is done in R using the [fgsea](https://bioconductor.org/packages/release/bioc/html/fgsea.html) package.
