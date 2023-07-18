# sangerseq2taxonomy

***Under condstruction - Use at your own risk***

This repo contains R scripts for processing Sanger sequencing <code>.ab1</code> files, and for generating descriptive statistic tables regarding the QC steps implemented.

Expected input: 
- 16S rRNA gene sequences in <code>.ab1</code> format
- Designed specifically for V3-V6 region amplicons sequenced 5'-3' backwards from V6 to V3 (i.e., reverse complement sequences)


## Quick start
### Step 1: 'abif_fasta2'
```r
#Source the following function from this repository
source("https://raw.githubusercontent.com/bdaisley/sangerseq2taxonomy/main/R_functions/functions-abif_fasta2.R")

#Copy  path of the folder where <code>.ab1</code> files are located.
#Example on a Windows-based system: "C:/bdaisley/sanger_files/2023_07_15"   # Ensure only forward slashes (/) in path

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Process and quality filter <code>.ab1</code> files of interest
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
abif_fasta2(folder="C:/bdaisley/sanger_files/2023_07_15",  
            export_html=TRUE,
            export_csv=TRUE,
            export_fasta=TRUE,
            export_fasta_revcomp=FALSE,
            verbose=TRUE,
            exclude=NULL,
            quality_cutoff = NULL,
            sliding_window_size = 15,
            date=NULL)

```

Inspect data via CSV files (containing pass/fail sequences) and HTML interactive [reactable](https://github.com/glin/reactable) output (see below)

<img src="https://github.com/bdaisley/sangerseq2taxonomy/blob/main/sangerseq2taxonomy.gif?raw=true" align="center" />


### Step 2: 'assign_taxonomy'
```r
#Source the following function from this repository
source("https://raw.githubusercontent.com/bdaisley/sangerseq2taxonomy/main/R_functions/functions-assign_taxonomy.R")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Assign taxonomy
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
assign_taxonomy(folder="C:/bdaisley/sanger_files/2023_07_15", #path should be same as in Step 1
                export_csv=TRUE,
                verbose=TRUE,
                skip_search=FALSE,
                quick_search=TRUE)

```

### Step 3: 'make_library'
```r
#Source the following function from this repository
source("https://raw.githubusercontent.com/bdaisley/sangerseq2taxonomy/main/R_functions/functions-make_library.R")

#Note: The input file to make a new library should be the PASS version of the 'assign_taxonomy' CSV output from Step 2
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Make new library file
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
make_library(new_lib_csv="C:/bdaisley/sanger_files/2023_07_15/output/assign_taxonomy_output_PASS___2023_07_15.csv",
             old_lib_csv=NULL,
             include_warnings=FALSE)
```

Inspect data via CSV files and HTML interactive [reactable](https://github.com/glin/reactable) output (see below)

<img src="https://github.com/bdaisley/sangerseq2taxonomy/blob/main/sangerseq2taxonomy_step3.gif?raw=true" align="center" />



More examples to come on usage of function...
