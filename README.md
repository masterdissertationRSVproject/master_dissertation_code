# master_dissertation_code
This repository contains the scripts used for analysing Celloseq data to map the RSV transcriptome at the single-cell level

All of the pre-processing step (Merge, internal adaptor QC, Porechop, and demultiplex) was run in CCB cluster, thus following the code available at https://github.com/MarioniLab/CELLOseq/ 

Below is the description of each script in this repository: 

Demultiplex True:False reads.R : this code was used to generate the histogram containing how many reads were ambigously mapped to each barcode during demultiplexing process.

FACS plot.R: this code was used to generate dot plots from .fcs data generated after FACS, displaying raw data collected for each event.

Gene Length vs Coverage for barcode 67.R : this code was used to examine the relationship between RSV gene length and coverage. This was also used to count the correlation coefficient for that data. The data was obtained from IGV after map all the reads from barcode 67 towards the RSV gene annotation file.

Plots for 4 conditions porechop.R : this code was used to generate various plot depicting read counts per cell, either mapped to the human or RSV genome. 

RSV Gene Expression barcode 67.R : this code was used to count the rsv gene expression for reads derived from barcode 67 after normalise the counts into CPM.

Z-Score Analysis : this code was used to perform statistical analysis including negative binomial model, Kolmogorov-Smirnov, Mann-Whitney U, and Z-score analysis for the reads derived from barcode 67
