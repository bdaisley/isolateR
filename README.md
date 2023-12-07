## isolateR: Automated processing of Sanger sequencing data, taxonomic profiling, and generation of microbial strain libraries

<p align="center"><img src="https://github.com/bdaisley/isolateR/blob/main/man/figures/isolateR_overview.jpg?raw=true" width="1500"></p>

***Beta version - Use at your own risk***

isolateR aims to enhance microbial isolation workflows and support the identification of novel taxa. It addresses the challenges of manual Sanger sequencing data processing and limitations of conventional BLAST searches, crucial for identifying microorganisms and creating strain libraries. The package offers a streamlined three-step process that automates quality trimming Sanger sequence files, taxonomic classification via global alignment against type strain databases, and efficient strain library creation based on customizable sequence similarity thresholds. It features interactive HTML output tables for easy data exploration and optional tools for generating phylogenetic trees to visualize microbial diversity.


- The expected input is Sanger sequence <code>.ab1</code> files containing taxonomic marker sequences.
- The pipeline is currently optimized for the following taxonomic markers:
  - <b>16S rRNA</b> (bacteria/archaea) 
  - <b>18S rRNA</b> (fungi)
  - <b>ITS region</b> (fungi)
  - <b>cpn60</b> (bacteria/archaea)

## Installation
```r
# install.packages("devtools")
devtools::install_github("bdaisley/isolateR")
```


## Quick start

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
<img src="https://github.com/bdaisley/isolateR/blob/main/man/figures/isolateR_step1_output.gif?raw=true" align="center" />


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
                    phylum_cutoff=75.0,
                    class_cutoff=78.5,
                    order_cutoff=82.0,
                    family_cutoff=86.5,
                    genus_cutoff=96.5,
                    species_cutoff=98.7)
```
```r			
# Parameters:
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# input			CSV file containing PASS sequences from isoQC step
# export_html		Toggle (TRUE/FALSE). Default=TRUE export results in CSV table.
# export_csv		Toggle (TRUE/FALSE). Default=TRUE export results in CSV table.
# db			Database for taxonomic classification ("16S","18S","ITS", or "cpn60")
# quick_search		Toggle (TRUE/FALSE) Default=FALSE performs comprehensive database search.
# phylum_cutoff		Similarity cutoff for Phylum rank demarcation (0-1)
# class_cutoff		Similarity cutoff for Class rank demarcation (0-1)
# order_cutoff		Similarity cutoff for Order rank demarcation (0-1)
# family_cutoff		Similarity cutoff for Family rank demarcation (0-1)
# genus_cutoff		Similarity cutoff for Genus rank demarcation (0-1)
# species_cutoff	Similarity cutoff for Species rank demarcation (0-1)

```



### Step 3: <code>isoLIB</code> - Generate strain library

This function creates a strain library by grouping closely related strains of interest based on sequence similarity.

For adding new sequences to an already-established strain library, specify the file path of the older strain library using the 'old_lib_csv" parameter.

- Note: The input file to make a new library should be the CSV output from 'isoTAX' in Step 2.


```r
#Specify location of CSV output from isoTAX in Step 2
fpath3 <- file.path(fpath1, "isolateR_output/02_isoTAX_results.csv")

isoLIB.S4 <- isoLIB(input=fpath3,
		    old_lib_csv=NULL,
		    include_warnings=FALSE,
		    strain_group_cutoff=0.995)
```
```r
# Parameters:
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# input				CSV file containing PASS sequences from isoTAX step
# old_lib_csv			If adding to existing library, provide 'isoLIB' output (.CSV extension) from past run. 
# include_warnings		Toggle (TRUE/FALSE) Set to TRUE to keep sequences with warnings from 'isoTAX' step.
# strain_group_cutoff		Similarity cutoff (0-1) for delineating strain groups. (1 = 100% identical/0.95=5.0% difference/etc.)
```


Inspect data via CSV files and HTML interactive [reactable](https://github.com/glin/reactable) output (see below)

<img src="https://github.com/bdaisley/isolateR/blob/main/man/figures/isolateR_step3_output.gif?raw=true" align="center" />


More examples on usage of functions to come...

