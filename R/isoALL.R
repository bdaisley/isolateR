#' @title Perform all commands in one step.
#' @name isoLIB
#' @rdname isoLIB
#' @description This function effectively wraps isoQC, isoTAX, and isoLIB steps into a single command for convenience. Input can be a single directory
#' or a list of directories to process at once. If multiple directories are provided, the resultant libraries can be sequentially merged together 
#' by toggling the parameter 'merge=TRUE'. All other respective parameters from the wrapped functions can be passed through this command.
#' . The The respective input parameters from
#' the wrappred  can be passed through this command with exception of the .creates a strain library by grouping closely related strains of interest based on sequence similarity.
#' For adding new sequences to an already-established strain library, specify the .CSV file path of the older strain library using the 'old_lib_csv" parameter.
#' @export
#' @param input Directory path(s) containing .ab1 files. If more than one, provivde as list (e.g. 'input=c("/path/to/directory1","/path/to/directory2")')
#' @param export_html (Default=TRUE) Output the results as an HTML file
#' @param export_csv (Default=TRUE) Output the results as a CSV file.
#' @param export_fasta (Default=TRUE) Output the sequences in a FASTA file.
#' @param export_fasta_revcomp (Default=FALSE) Output the sequences in reverse complement form in a fasta file. This is useful in cases where sequencing was done using the reverse primer and thus the orientation of input sequences needs reversing.
#' @param verbose (Default=FALSE) Output progress while script is running.
#' @param files_manual (Default=NULL) For testing purposes only. Specify a list of files to run  as filenames without extensions, rather than the whole directory format. Primarily used for testing, use at your own risk.
#' @param exclude (Default=NULL) For testing purposes only. Excludes files of interest from input directory.
#' @param min_phred_score (Default=20) Do not accept trimmed sequences with a mean Phred score below this cutoff
#' @param min_length (Default=200) Do not accept trimmed sequences with sequence length below this number
#' @param sliding_window_cutoff (Default=NULL) Quality trimming parameter (M2) for wrapping SangerRead function in sangeranalyseR package. If NULL, implements auto cutoff for Phred score (recommended), otherwise set between 1-60.
#' @param sliding_window_size (Default=15) Quality trimming parameter (M2) for wrapping SangerRead function in sangeranalyseR package. Recommended range between 5-30.
#' @param date Set date "YYYY_MM_DD" format. If NULL, attempts to parse date from .ab1 file
#' @param quick_search (Default=FALSE) Whether or not to perform a comprehensive database search (i.e. optimal global alignment).
#' If TRUE, performs quick search equivalent to setting VSEARCH parameters "--maxaccepts 100 --maxrejects 100".
#' If FALSE, performs comprehensive search equivalent to setting VSEARCH parameters "--maxaccepts 0 --maxrejects 0"
#' @param db (Default="16S") Select database option(s) including "16S" (for searching against the NCBI Refseq targeted loci 16S rRNA database),
#' "ITS" (for searching against the NCBI Refseq targeted loci ITS  database. For combined databases in cases where input sequences are dervied from
#' bacteria and fungi, select "16S|ITS".
#' @param iddef Set pairwise identity definition as per VSEARCH definitions (Default=2, and is recommended for highest taxonomic accuracy)
#' (0) CD-HIT definition: (matching columns) / (shortest sequence length).
#' (1) Edit distance: (matching columns) / (alignment length).
#' (2) Edit distance excluding terminal gaps (default definition).
#' (3) Marine Biological Lab definition counting each gap opening (internal or terminal) as a single mismatch, whether or not the gap was extended: 1.0- ((mismatches + gap openings)/(longest sequence length)).
#' (4) BLAST definition, equivalent to --iddef 1 for global pairwise alignments.
#' @param phylum_cutoff Percent cutoff for phylum rank demarcation
#' @param class_cutoff Percent cutoff for class rank demarcation
#' @param order_cutoff Percent cutoff for order rank demarcation
#' @param family_cutoff Percent cutoff for family rank demarcation
#' @param genus_cutoff Percent cutoff for genus rank demarcation
#' @param species_cutoff Percent cutoff for species rank demarcation
#' @param include_warnings (Default=FALSE) Whether or not to keep sequences with poor alignment warnings from Step 2 'isoTAX' function. Set TRUE to keep warning sequences, and FALSE to remove warning sequences.
#' @param strain_group_cutoff (Default=0.995) Similarity cutoff (0-1) for delineating between strain groups. 1 = 100% identical/0.995=0.5% difference/0.95=5.0% difference/etc.
#' @seealso \code{\link{isoQC}}, \code{\link{isoTAX}}, \code{\link{isoLIB}}
#' @return Returns a list of \code{\link{class-isoLIB}} class objects.
#' @examples
#' #Set path to directory containing example .ab1 files
#' fpath1 <- system.file("extdata/abif_examples/rocket_salad", package = "isolateR")
#'
#' #Run isoALL function with default settings
#' isoALL(input=fpath1)

isoALL <- function(input=NULL,
                   export_html=TRUE,
                   export_csv=TRUE,
                   export_fasta=TRUE,
                   export_fasta_revcomp=FALSE,
                   quick_search=FALSE,
                   db="16S",
                   iddef=2,
                   phylum_cutoff=75.0,
                   class_cutoff=78.5,
                   order_cutoff=82.0,
                   family_cutoff=86.5,
                   genus_cutoff=96.5,
                   species_cutoff=98.7,
                   include_warnings=FALSE,
                   strain_group_cutoff = 0.995,
                   merge=FALSE){
  
  #Initial input checks
  #--------------------------------------------------
  if(is.null(input)) stop("Name of folder not supplied. Please enter the path to the folder containing .ab1 sequence files.")
  if(length(input) < 2 & merge==TRUE) stop("Only 1 input directory provided, cannot merge files. Please enter a list of directories, or set 'merge=FALSE' to proceed.")
  export_csv <- TRUE
  #--------------------------------------------------
  
  isoALL.list <- list()
  
  exists.arg <- file.exists(paste(input[1], "/isolateR_output/03_isoLIB_results_merged.csv", sep=""))
  
  if((exists.arg & merge==TRUE)){
    message(cat(paste0("\n", "\033[97;", 40, "m", 
                       "Previously merged library found in first directory. Starting at second directory and continuing to build.", 
                       "\033[0m")))
    for(i in 2:length(input)){
      #::::::::::::::::::::
      # Step 1) isoQC
      #::::::::::::::::::::
      isoQC.S4 <- isoQC(input=input[i],
                        export_html=export_html,
                        export_csv=TRUE,
                        export_fasta=export_fasta,
                        export_fasta_revcomp=export_fasta_revcomp,
                        verbose=FALSE,
                        exclude=NULL,
                        min_phred_score=20,
                        min_length=200,
                        sliding_window_cutoff=NULL,
                        sliding_window_size = 15,
                        date=NULL,
                        files_manual=NULL)
      
      #::::::::::::::::::::
      # Step 2) isoTAX
      #::::::::::::::::::::
      isoTAX.S4 <- isoTAX(input=paste(input[i], "/isolateR_output/01_isoQC_trimmed_sequences_PASS.csv", sep=""),
                          export_html=export_html,
                          export_csv=TRUE,
                          quick_search=quick_search,
                          db=db,
                          iddef=iddef,
                          phylum_cutoff=phylum_cutoff,
                          class_cutoff=class_cutoff,
                          order_cutoff=order_cutoff,
                          family_cutoff=family_cutoff,
                          genus_cutoff=genus_cutoff,
                          species_cutoff=species_cutoff)
      
      #::::::::::::::::::::
      # Step 3) isoLIB
      #::::::::::::::::::::
      
      old_lib_csv.x <- paste(input[(i-1)], "/isolateR_output/03_isoLIB_results_merged.csv", sep="")
      
      isoLIB.S4 <- isoLIB(input=paste(input[i], "/isolateR_output/02_isoTAX_results.csv", sep=""),
                          old_lib_csv=old_lib_csv.x,
                          export_html=export_html,
                          export_csv=TRUE,
                          include_warnings=include_warnings,
                          strain_group_cutoff = strain_group_cutoff,
                          phylum_cutoff=phylum_cutoff,
                          class_cutoff=class_cutoff,
                          order_cutoff=order_cutoff,
                          family_cutoff=family_cutoff,
                          genus_cutoff=genus_cutoff,
                          species_cutoff=species_cutoff)
      
      isoALL.list <- c(isoALL.list, isoLIB.S4)
    }
  }
  
  if(!(exists.arg & merge==TRUE)){
    for(i in 1:length(input)){
      #::::::::::::::::::::
      # Step 1) isoQC
      #::::::::::::::::::::
      isoQC.S4 <- isoQC(input=input[i],
                        export_html=export_html,
                        export_csv=TRUE,
                        export_fasta=export_fasta,
                        export_fasta_revcomp=export_fasta_revcomp,
                        verbose=FALSE,
                        exclude=NULL,
                        min_phred_score=20,
                        min_length=200,
                        sliding_window_cutoff=NULL,
                        sliding_window_size = 15,
                        date=NULL,
                        files_manual=NULL)
      
      #::::::::::::::::::::
      # Step 2) isoTAX
      #::::::::::::::::::::
      isoTAX.S4 <- isoTAX(input=paste(input[i], "/isolateR_output/01_isoQC_trimmed_sequences_PASS.csv", sep=""),
                          export_html=export_html,
                          export_csv=TRUE,
                          quick_search=quick_search,
                          db=db,
                          iddef=iddef,
                          phylum_cutoff=phylum_cutoff,
                          class_cutoff=class_cutoff,
                          order_cutoff=order_cutoff,
                          family_cutoff=family_cutoff,
                          genus_cutoff=genus_cutoff,
                          species_cutoff=species_cutoff)
      
      #::::::::::::::::::::
      # Step 3) isoLIB
      #::::::::::::::::::::
      old_lib_csv.x <- NULL
      
      if(merge==TRUE & i > 1){old_lib_csv.x <- paste(input[(i-1)], "/isolateR_output/03_isoLIB_results.csv", sep="")}
      if(merge==TRUE & i > 2){old_lib_csv.x <- paste(input[(i-1)], "/isolateR_output/03_isoLIB_results_merged.csv", sep="")}
      
      isoLIB.S4 <- isoLIB(input=paste(input[i], "/isolateR_output/02_isoTAX_results.csv", sep=""),
                          old_lib_csv=old_lib_csv.x,
                          export_html=export_html,
                          export_csv=TRUE,
                          include_warnings=include_warnings,
                          strain_group_cutoff = strain_group_cutoff,
                          phylum_cutoff=phylum_cutoff,
                          class_cutoff=class_cutoff,
                          order_cutoff=order_cutoff,
                          family_cutoff=family_cutoff,
                          genus_cutoff=genus_cutoff,
                          species_cutoff=species_cutoff)
      
      isoALL.list <- c(isoALL.list, isoLIB.S4)
    }
  }
}
