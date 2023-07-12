# sangerseq2taxonomy

***Under condstruction - Use at your own risk***

This repo contains R scripts for processing Sanger sequencing <code>.ab1</code> files, and for generating descriptive statistic tables regarding the QC steps implemented.

Expected input: 
- 16S rRNA gene sequences in <code>.ab1</code> format
- Designed specifically for V3-V6 region amplicons sequenced 5'-3' backwards from V6 to V3 (i.e., reverse complement sequences)


## Quick start
Step 1: Source the following function from this repository
```r
source("https://raw.githubusercontent.com/bdaisley/sangerseq2taxonomy/main/R_functions/functions-abif_fasta2.R")
```
Step 2: Set working paths to where .ab1 files are located
```r
folder.name <- "2023_07_06" #Note: this is the folder name only, not entired path
setwd(paste("~/Sanger_sequencing_results/", folder.name, sep="")) #Note: this the entire path *before* the folder name above
sanger.path <- getwd()
```
Step 3: Run function of interest to process .ab1 files
```r
abif_fasta2(folder=sanger.path, 
            reversecomp=TRUE,
            export_html=TRUE,
            export_fasta=TRUE,
            verbose=TRUE,
            output="V3-V6seq.fasta")
```
Step 4: Inspect data via interactive reactable output

<img src="https://github.com/bdaisley/sangerseq2taxonomy/blob/main/sangerseq2taxonomy.gif?raw=true" align="center" />

More to come on automated classification of taxonomy...