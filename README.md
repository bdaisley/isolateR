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
Step 2: Set working paths to where <code>.ab1</code> files are located. Example on a Windows-based system:
```r
sanger.path <- "C:/Users/Brendan/Documents/Sanger_sequencing_results/2023_07_06"
```
Step 3: Run function of interest to process <code>.ab1</code> files
```r
abif_fasta2(folder=sanger.path,
            export_html=TRUE,
            export_csv=TRUE,
            export_fasta=TRUE,
            export_fasta_revcomp=TRUE,
			quality_cutoff = 25,
            sliding_window_size = 15
            verbose=TRUE)
```
Step 4: Inspect data via interactive reactable output

<img src="https://github.com/bdaisley/sangerseq2taxonomy/blob/main/sangerseq2taxonomy.gif?raw=true" align="center" />

More to come on automated classification of taxonomy...