# Cell_Spec_EPI

This repository houses the codes utilzied in the paper "Network Analysis of Enhancer-Promoter Interactions Highlights Cell-Type Specific Mechanisms of Transcriptional Regulation Variation"

This repository is currently still in the process of uploading the codes and documentation used in the study. 

If you have any questions regarding the codes used please reach out to jhk148@dls.rutgers.edu

### ATCRCP Scripts
These 2 files, 1 for paired end data and 1 for single ended data, are utilized to convert the raw ATACseq, CUT&RUN, or ChIPseq fastq files into processed files for use in generating ABC model predicted enhancer promoter interactions. The files process 1 assay type and 1 replicate at a time.

### ABC_CommandSetup
This python script takes in all the above processed data and generates the bash code that is to be run for generating ABC model enhancer promoter interactions



