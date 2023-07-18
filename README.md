## isolateR: sangerseq to taxonomy pipeline for automating microbial isolation workflows and the seamless generation of novel strain libraries

<p align="center"><img src="https://github.com/bdaisley/isolateR/blob/main/isolateR_overview.jpg?raw=true" width="1500"></p>

***Under construction - Use at your own risk***

This repo contains R scripts for processing Sanger sequencing <code>.ab1</code> files, and for generating descriptive statistic tables regarding the QC steps implemented.

Expected input: 
- 16S rRNA gene sequences in <code>.ab1</code> format
- Designed specifically for V3-V6 region amplicons sequenced 5'-3' backwards from V6 to V3 (i.e., reverse complement sequences)


## Quick start
### Step 1: <code>abif_fasta2</code>
```r
#Source the following function from this repository
source("https://raw.githubusercontent.com/bdaisley/isolateR/main/R_functions/functions-abif_fasta2.R")

#Copy  path of the folder where <code>.ab1</code> files are located.
#Example on a Windows-based system: "C:/bdaisley/sanger_files/2023_07_06"   # Ensure only forward slashes (/) in path

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Process and quality filter <code>.ab1</code> files of interest
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
abif_fasta2(folder="C:/bdaisley/sanger_files/2023_07_06",  # Folder containing .ab1 files
            export_html=TRUE,                              # Interactive display of output
            export_csv=TRUE,                               # Exports separate files for PASS/FAIL seqs, needed for next step 'assign_taxonomy'
            export_fasta=TRUE,                             # Exports separate files for PASS/FAIL seqs
            export_fasta_revcomp=FALSE,                    # Exports reverse complement of fasta sequences
            verbose=TRUE,                                  # Toggle messages in R console on/off
            quality_cutoff = NULL,                         # NULL by default implements auto cutoff (recommended)
            sliding_window_size = 15,                      # For quality trimming steps. Default = 15 (recommended)
            date=NULL)                                     # Set date "YY_MM_DD" format. If NULL, attempts to parse date from .ab1 file

```

Inspect data via CSV files (containing pass/fail sequences) and HTML interactive [reactable](https://github.com/glin/reactable) output (see below)
isolateR_step1_output.gif
<img src="https://github.com/bdaisley/isolateR/blob/main/isolateR_step1_output.gif?raw=true" align="center" />


### Step 2: <code>assign_taxonomy</code>
```r
#Source the following function from this repository
source("https://raw.githubusercontent.com/bdaisley/isolateR/main/R_functions/functions-assign_taxonomy.R")

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Assign taxonomy
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
assign_taxonomy(folder="C:/bdaisley/sanger_files/2023_07_06",  # Path should be same as in Step 1 above
                export_csv=TRUE,                               # Exports taxonomic assignment CSV file, needed for input of next step 'make_library'
                verbose=TRUE,                                  # Toggle messages in R console on/off
                skip_search=FALSE,                             # FALSE by default. If TRUE, skips database searching step (requires presence of output files from previous run)
                quick_search=TRUE,                             # TRUE for convenience/speed in this example. FALSE by default, which performs comprehensive database searching (recommended for confident taxonomy) 
                add_fungi_db=FALSE)                            # Set TRUE if any .ab1 files contain ITS sequences. FALSE by default, which only searches against NCBI bacterial 16S rRNA database

```

### Step 3: <code>make_library</code>
```r
#Source the following function from this repository
source("https://raw.githubusercontent.com/bdaisley/isolateR/main/R_functions/functions-make_library.R")

#Note: The input file to make a new library should be the PASS version of the 'assign_taxonomy' CSV output from Step 2
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#Make new library file
#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
make_library(new_lib_csv="C:/bdaisley/sanger_files/2023_07_06/output/assign_taxonomy_output_PASS___2023_07_06.csv",  # Output .csv file from 'assign_taxonomy' in Step 2 above. Note: input is a file, not a folder.
             old_lib_csv=NULL,                                                                                       # If adding to existing library, provide a previously generated 'make_library' .csv output file. 
             include_warnings=FALSE)                                                                                 # Set to TRUE to keep sequences with poor alignment warnings from 'assign_taxonomy' in Step 2 above. FALSE by default.
```

Inspect data via CSV files and HTML interactive [reactable](https://github.com/glin/reactable) output (see below)

<img src="https://github.com/bdaisley/isolateR/blob/main/isolateR_step3_output.gif?raw=true" align="center" />



More examples on usage of functions to come...

