# Scripts used as part of the pipeline for the analysis of mitochondrial reads
We run this pipeline using a combination of in-house shell, python, and R scripts.

## Shell scripts

These are use as part of the main pipeline, controlled by the two master shell scripts (master.sh and master.pairedend.sh, depending on the nature of the input data: either single-end or paired-end). All other shell scripts are called by one of the two master scripts, apart from accuracy.sh that is used on simulated datasets to quantify aligners accuracy and collect other metrics.

## Python scripts

These are used when the handling of large datasets is required.

## R script

This script is used for plotting data using the [karyoplotR](https://bioconductor.org/packages/release/bioc/html/karyoploteR.html) package.
