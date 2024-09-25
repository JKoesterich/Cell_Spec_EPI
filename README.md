# Cell_Spec_EPI

This repository houses the codes utilzied in the paper "Network Analysis of Enhancer-Promoter Interactions Highlights Cell-Type Specific Mechanisms of Transcriptional Regulation Variation"

This repository is currently still in the process of uploading the codes and documentation used in the study. 

If you have any questions regarding the codes used please reach out to jhk148@dls.rutgers.edu

### ATCRCP Scripts
These 2 script files, 1 for paired end data and 1 for single ended data, are utilized to convert the raw ATACseq, CUT&RUN, or ChIPseq fastq files into processed files for use in generating ABC model predicted enhancer promoter interactions. The files process 1 assay type and 1 replicate at a time.

### ABC_CommandSetup Script
The python script with this name takes in all the above processed data and generates the bash code that is to be run for generating ABC model enhancer promoter interactions

## Generating Cluster names
After the ABC model is run using the above generated commands, we take the EnhancerPredictionsFull.txt output file and filter out any interactions below 0.02. Then we remove any lines where the target gene ID has less than 1 TPM in that cell type's processed RNA file (provided from original studies paper). We then select only the columns for the bed position of the enhancer, the gene ID of the promoter, and the interaction score from the resulting file and add the cell type label to the end of the line. This is done by using awk to print the line with the added cell type label text as a new 5th column

### Adding Sub-Structure label of EPI regulation
Using the 6C_Concat_CNum_PerClusterIdentifyer python script, we add the sub-structure labels to the rows of EPIs. We first require that the file is presorted by the cell type label to make sure we can categorize the entire file in 1 read through and only have the sub-structure determined by the EPIs within the same cell type. This uses the frist 3 columns to identify if the same enhancer location is found multiple time and uses the 4th column to determine if the gene name is seen multiple times. Then after a full cell type cluster is analyzed (either by observing the next entry having a different cell type cluster label or reaching the end of the file. The script goes through the stored EPIs and prints the full row along with concatenating the sub-structure label as an additional column at the end.

