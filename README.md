<p align="center"><img src="https://github.com/bdaisley/isolateR/blob/main/man/figures/isolateR.png?raw=true"></p>

## isolateR: Automated processing of Sanger sequencing data, taxonomic profiling, and generation of microbial strain libraries

<p align="center"><img src="https://github.com/bdaisley/isolateR/blob/main/man/figures/isolateR_overview1.jpg?raw=true" width="1500"></p>

<b>Update July 2024:</b> isolateR is now published in [Bioinformatics](https://academic.oup.com/bioinformatics/article/40/7/btae448/7712426)!

isolateR aims to enhance microbial isolation workflows and support the identification of novel taxa. It addresses the challenges of manual Sanger sequencing data processing and limitations of conventional BLAST searches, crucial for identifying microorganisms and creating strain libraries. The package offers a streamlined three-step process that automates quality trimming Sanger sequence files, taxonomic classification via global alignment against type strain databases, and efficient strain library creation based on customizable sequence similarity thresholds. It features interactive HTML output tables for easy data exploration and optional tools for generating phylogenetic trees to visualize microbial diversity. 


- The expected input is Sanger sequence <code>.ab1</code> files containing taxonomic marker sequences.
- The pipeline is currently optimized for the following taxonomic markers:
  - <b>16S rRNA</b> (bacteria/archaea) 
  - <b>18S rRNA</b> (fungi)
  - <b>ITS region</b> (fungi)
  - <b>cpn60</b> (bacteria/archaea)

## Installation

### Install via GitHub (Development version)
```r
#Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#Install the required Bioconductor dependencies
BiocManager::install(c("Biostrings", "msa", "sangerseqR"), update=FALSE)

#Install isolateR
remotes::install_github("bdaisley/isolateR")
```
### Install via Conda
```shell
conda create --name isolateR r-isolater=1.0.1-0 -c bdaisley -c conda-forge -c bioconda
```

## Quick start
The one-step command <code>isoALL</code> wraps the three main functions below to quickly process .ab1 files in batch. Simply specify the folder(s) containing .ab1 sequence files you want to process.

#### Example with a single input folder
```r
library(isolateR)
isoALL.S4 <- isoALL(input="/path/to/folder_containing_ab1_files")
```
#### Example with multiple input folders + merging results
```r
library(isolateR)
folder_list <- c("/path/to/folder_containing_ab1_files1",
                 "/path/to/folder_containing_ab1_files2",
                 "/path/to/folder_containing_ab1_files3")

isoALL.S4 <- isoALL(input=folder_list, merge=TRUE)
```

## Overview of the 3 main functions

### Step 1: <code>isoQC</code> - Automated quality trimming of sequences

This function loads in ABIF files (.ab1 extension) and performs automatic quality trimming in batch mode.

Reminder on Windows-based operating systems: Ensure only forward slashes (/) used in your path when setting directory.

```r
library(isolateR)

#Set path of directory where the .ab1 files. In this case, using example dataset in R
fpath1 <- system.file("extdata/abif_examples/rocket_salad", package = "isolateR")

isoQC.S4 <- isoQC(input=fpath1,
                  export_html=TRUE,
                  export_csv=TRUE,
                  export_fasta=TRUE,
                  verbose=FALSE,
                  min_phred_score = 20,
                  min_length = 200,
                  sliding_window_cutoff = NULL,
                  sliding_window_size = 15,
                  date=NULL)
```
```r
# Parameters:
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# input				Path of directory containing input .ab1 files
# export_html			Toggle export of results in interactive HTML table.(TRUE/FALSE).
# export_csv			Toggle export of PASS/FAIL sequence results in CSV format (TRUE/FALSE).
# export_fasta			Toggle export of PASS/FAIL sequences in FASTA format(TRUE/FALSE).
# verbose			Toggle checkpoint messages in R console (TRUE/FALSE).
# min_phred_score		Do not accept trimmed seqs with phred score cutoff below this number. (Default=20)
# min_length			Do not accept trimmed seqs with sequence length below this number
# sliding_window_cutoff		For quality trimming steps. NULL by default implements auto cutoff (recommended).
# sliding_window_size		For quality trimming steps. (Default= 15)
# date				Set date "YYYY_MM_DD" format. (Default=NULL) attempts to parse date from .ab1 file.
```

The exported CSV files containing PASS/FAIL sequences based on quality thresholds include:
- "01_isoQC_trimmed_sequences_PASS.csv"
- "01_isoQC_trimmed_sequences_FAIL.csv"

Descriptive statistics regarding QC steps implemented can be inspected via interactive HTML tables in the [reactable](https://github.com/glin/reactable) output (see below)
isolateR_step1_output.gif
<img src="https://github.com/bdaisley/isolateR/blob/main/man/figures/isolateR_step1_output1.gif?raw=true" align="center" />


### Step 2: <code>isoTAX</code> - Assign taxonomy

This function performs taxonomic classification by searching query Sanger sequences against specified database of interest. Takes CSV input files, extracts FASTA-formatted query sequences and performs global alignment against specified database of interest via Needleman-Wunsch algorithm by wrapping the --usearch_global command implemented in VSEARCH. Default taxonomic rank cutoffs for 16S rRNA gene sequences are based on Yarza et al. 2014, Nat Rev Microbiol.

- The input for this step is expected to be the .CSV file exported in the previous step (e.g. "01_isoQC_trimmed_sequences_PASS.csv")

- Note: It is possible for users to manually add back in failed sequences by appending rows of interest from the fail .CSV output to the pass .CSV or by combining them in a new .CSV document altogether. In such a case, the column names and dimensions must be identical to the original output.

```r
#Specify location of CSV output from 'isoQC' step containing quality trimmed sequences
fpath2 <- file.path(fpath1, "isolateR_output/01_isoQC_trimmed_sequences_PASS.csv")

isoTAX.S4 <- isoTAX(input=fpath2,
                    export_html=TRUE,
                    export_csv=TRUE,
                    db="16S",
                    quick_search=TRUE,
                    phylum_threshold=75.0,
                    class_threshold=78.5,
                    order_threshold=82.0,
                    family_threshold=86.5,
                    genus_threshold=96.5,
                    species_threshold=98.7)
```
```r			
# Parameters:
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# input			CSV file containing PASS sequences from isoQC step
# export_html		Toggle (TRUE/FALSE). Default=TRUE export results in CSV table.
# export_csv		Toggle (TRUE/FALSE). Default=TRUE export results in CSV table.
# db			Database for taxonomic classification ("16S","18S","ITS", or "cpn60")
# quick_search		Toggle (TRUE/FALSE) Default=FALSE performs comprehensive database search.
# phylum_threshold	Similarity threshold for Phylum rank demarcation (0-100)
# class_threshold	Similarity threshold for Class rank demarcation (0-100)
# order_threshold	Similarity threshold for Order rank demarcation (0-100)
# family_threshold	Similarity threshold for Family rank demarcation (0-100)
# genus_threshold	Similarity threshold for Genus rank demarcation (0-100)
# species_threshold	Similarity threshold for Species rank demarcation (0-100)

```

<img src="https://github.com/bdaisley/isolateR/blob/main/man/figures/isolateR_step2_output1.gif?raw=true" align="center" />

### Step 3: <code>isoLIB</code> - Generate strain library

This function creates a strain library by grouping closely related strains of interest based on sequence similarity.

For adding new sequences to an already-established strain library, specify the file path of the older strain library using the 'old_lib_csv" parameter.

- Note: The input file to make a new library should be the CSV output from 'isoTAX' in Step 2.


```r
#Specify location of CSV output from isoTAX in Step 2
fpath3 <- file.path(fpath1, "isolateR_output/02_isoTAX_results.csv")

isoLIB.S4 <- isoLIB(input=fpath3,
		    old_lib_csv=NULL,
		    group_cutoff=0.995,
                    include_warnings=FALSE)
```
```r
# Parameters:
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# input				CSV file containing PASS sequences from isoTAX step
# old_lib_csv			If adding to existing library, provide 'isoLIB' output (.CSV extension) from past run. 
# group_cutoff			Similarity cutoff (0-1) for delineating strain groups. (1 = 100% identical/0.95=5.0% difference/etc.)
# include_warnings		Toggle (TRUE/FALSE) Set to TRUE to keep sequences with warnings from 'isoTAX' step.
```


Inspect data via CSV files and HTML interactive [reactable](https://github.com/glin/reactable) output (see below)

<img src="https://github.com/bdaisley/isolateR/blob/main/man/figures/isolateR_step3_output1.gif?raw=true" align="center" />


More examples on usage of functions to come...


## Citation

Daisley B., Vancuren S.J., Brettingham D.J.L., Wilde J., Renwick S., Macpherson C., Good D.A., Botschner A.J., Yen S., Hill J.E., Sorbara M.T., Allen-Vercone E. (2024). isolateR: an R package for generating microbial libraries from Sanger sequencing data. Bioinformatics 40(7):btae448. (https://doi.org/10.1093/bioinformatics/btae448)


