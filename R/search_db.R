#' @name search_db
#' @title Perform global alignment pairwise identity search using VSEARCH and type strain database of interest.
#' @description Performs global alignment between FASTA-formatted query sequences and the specified database of interest.
#' Uses the Needleman-Wunsch algorithm by wrapping the --usearch_global command implemented in VSEARCH.
#' @export
#' @param query.path Path of FASTA-formatted query sequence file.
#' @param uc.out Path of UC-formatted results output table.
#' @param b6.out Path of blast6-formatted results output table.
#' @param path Working path directory (Default is set to current working directory via 'getwd()'
#' @param quick_search (Default=FALSE) Whether or not to perform a comprehensive database search (i.e. optimal global alignment).
#' If TRUE, performs quick search equivalent to setting VSEARCH parameters "--maxaccepts 100 --maxrejects 100".
#' If FALSE, performs comprehensive search equivalent to setting VSEARCH parameters "--maxaccepts 0 --maxrejects 0"
#' Note: This option is provided for convenience and rough approximation of taxonomy only, set to FALSE for accurate % pairwise identity results.
#' @param db Optional: Select any of the standard database option(s) including "16S" (for searching against the NCBI Refseq targeted loci 16S rRNA database),
#' "ITS" (for searching against the NCBI Refseq targeted loci ITS  database. For combined databases in cases where input sequences are dervied from
#' bacteria and fungi, select "16S|ITS". Setting to anything other than db=NULL or db="custom" causes 'db.path' parameter to be ignored.
#' @param db.path Path of FASTA-formatted database sequence file. Ignored if 'db' parameter is set to anything other than "custom"
#' @param vsearch_path Path of VSEARCH software if manually downloaded in a custom directory. If NULL (Default), will attempt automatic download.
#' @param keep_temp_files Toggle (TRUE/FALSE). If TRUE, temporary .uc and .b6o output files are kept from VSEARCH --uc and --blast6out commands, respectively. If FALSE, temporary files are removed.
#' @param iddef Set pairwise identity definition as per VSEARCH definitions (Default=2, and is recommended for highest taxonomic accuracy)
#' (0) CD-HIT definition: (matching columns) / (shortest sequence length).
#' (1) Edit distance: (matching columns) / (alignment length).
#' (2) Edit distance excluding terminal gaps (default definition).
#' (3) Marine Biological Lab definition counting each gap opening (internal or terminal) as a single mismatch, whether or not the gap was extended: 1.0- ((mismatches + gap openings)/(longest sequence length)).
#' (4) BLAST definition, equivalent to --iddef 1 for global pairwise alignments.
#' @return Returns a dataframe matching the UC-formatted output table from VSEARCH. Query sequences are automatically added to the final column.
#' Summary of column information. See VSEARCH documentation for more details.
#'- V1 = Record type of hit (H) or no hit (N)
#'- V2 = Ordinal number of the target sequence (based on input order, starting from zero). Set to '*' for N.
#'- V3 = Sequence length. Set to '*' for N.
#'- V4 = Percentage of similarity with the target sequence. Set to '*' for N.
#'- V5 = Match orientation + or -. . Set to '.' for N.
#'- V6 = Not used, always set to zero for H, or '*' for N.
#'- V7 = Not used, always set to zero for H, or '*' for N.
#'- V8 = Compact representation of the pairwise alignment using the CIGAR format (Compact Idiosyncratic Gapped Alignment Report): M (match/mismatch), D (deletion) and I (insertion). The equal sign '=' indicates that the query is identical to the centroid sequence. Set to '*' for N.
#'- V9 = Label of the query sequence. Equivalent to 'filename' slot of isolateR class objects (e.g. isoQC, isoTAX, isoLIB).
#'- V10 = Label of the target centroid sequence. Set to '*' for N.
#' @import dplyr
#' @importFrom seqinr write.fasta
#' @importFrom seqinr read.fasta
#' @importFrom R.utils gunzip
#' @importFrom utils download.file
#' @importFrom stringr str_subset
#' @seealso \code{\link{isoTAX}}
#' @examples
#' #Set path to directory containing example .ab1 files
#' fpath1 <- system.file("extdata/abif_examples/rocket_salad", package = "isolateR")
#'
#' #Run isoQC function with default settings
#' isoQC.S4 <- isoQC(input=fpath1)
#'
#' #Set path of CSV output file containing PASS sequences from isoQC step
#' fasta.path <- "01_isoQC_trimmed_sequences_PASS.fasta"
#'
#' #Set paths
#' output.path <- file.path(fpath1, "isolateR_output")
#'
#' #Run search_db function
#' uc.df <- search_db(query.path=fasta.path, path=output.path, quick_search=TRUE, db="16S")
#'
#' #Inspect results
#' uc.df[1:10,1:10]

search_db <- function(query.path = NULL,
                      uc.out = "VSEARCH_output.uc",
                      b6.out = "VSEARCH_output.b6o",
                      path = getwd(),
                      quick_search = FALSE,
                      db = NULL,
                      db_path = NULL,
                      vsearch_path=NULL,
                      keep_temp_files=FALSE,
                      iddef=2){

  #:::::::::::::::::::::::::::::::::
  #Paths and parameter checks
  #:::::::::::::::::::::::::::::::::
  db.path <- db_path
  if(is.null(db) & !is.null(db.path)){message('Setting db="custom" as `db_path` was specified without specifiying `db`')
    db <- "custom"
    }
  if(is.null(db) & is.null(db.path)) stop("...Aborting step....Both 'db.path' and 'db' cannot be set to NULL.
  Select any of the standard database option(s) including '16S' (for searching
  against the NCBI Refseq targeted loci 16S rRNA database), 'ITS' (for searching
  against the NCBI Refseq targeted loci ITS  database. For combined databases in
  cases where input sequences are dervied from bacteria and fungi, select '16S|ITS'.
  Setting to anything other than NULL causes 'db.path' parameter to be ignored.", call.=FALSE)

  #set paths--------------------------------------------------------------------------------------------
  setwd(path)

  #:::::::::::::::::::::::::::::::::
  #Download reference databases
  #:::::::::::::::::::::::::::::::::
  
  if(db!="custom"){
    db.splits <- unlist(stringr::str_split(db, pattern="\\|"))
    
    #Single database requested------------------------------------------------------
    if(length(db.splits) == 1){db.path <- get_db(db=db)}
    
    #Combine databases if requested-------------------------------------------------
    if(length(db.splits) > 1){
      db.splits.list <- lapply(db.splits, function(x){
        db.tmp <- Biostrings::readBStringSet(get_db(x))
        as.data.frame(cbind("names" = names(db.tmp), "seqs" = paste(db.tmp)))
      })
      db.splits.list.combined <- dplyr::bind_rows(db.splits.list)
      fasta.combined <- DNAStringSet(setNames(db.splits.list.combined$seqs, db.splits.list.combined$names))
      db.fasta.path <- paste(file.path(system.file("", package="isolateR"), "databases"), "/", paste(db.splits, collapse="_"), ".fna", sep="")
      #Write combined database to FASTA file
      Biostrings::writeXStringSet(fasta.combined, file=db.fasta.path)
      db.path <- db.fasta.path
    }
  }

  if(db=="custom"){
    if(is.null(db.path)){ stop("The database file specified does not exist: db.path=", db.path, ". Please provide the correct path to your database file.",call.=FALSE)}
    if(!file.exists(db.path)){ stop("The database file specified does not exist: db.path=", db.path, ". Please provide the correct path to your database file.",call.=FALSE)}
    db.path <- db.path
  }
    
  #:::::::::::::::::::::::::::
  #Download VSEARCH software
  #:::::::::::::::::::::::::::
  
  if(is.null(vsearch_path)){
    vsearch.path <- get_vsearch()
  } else {
    vsearch.path <- vsearch_path
  }
  
  #:::::::::::::::::::::
  # Search function
  #:::::::::::::::::::::
  message(cat(paste0("\033[97;", 40, "m","Searching query sequences against NCBI database", "\033[0m")))
  message(cat(paste0("\033[0;", 32, "m","This may take several minutes...", "\033[0m")))

  suppressWarnings(dir.create(file.path(path, "temp_vsearch")))
  path <- file.path(path, "temp_vsearch")
  uc.out <- file.path("temp_vsearch", uc.out)
  b6.out <- file.path("temp_vsearch", b6.out)

  if(quick_search==TRUE){
    #Set compatible database path for VSEARCH
    invisible(file.copy(db.path, path, overwrite = TRUE))
    db.path.x <- unlist(strsplit(db.path, '/'))[length(unlist(strsplit(db.path, '/')))]
    db.path.x <- file.path("temp_vsearch", db.path.x)
    system2(vsearch.path, paste(" --usearch_global ", query.path, " --db ", db.path.x, " --blast6out ", b6.out, " --uc ", uc.out, " --iddef ",iddef, " --id 0.7 -query_cov 0.95 --maxaccepts 100 --maxrejects 100 --top_hits_only --strand both", sep=""), stdout="", stderr="")
    #Clean up
    unlink(file.path(path, db.path.x),recursive=TRUE)
    }
  if(quick_search==FALSE){
    #Set compatible database path for VSEARCH
    invisible(file.copy(db.path, path, overwrite = TRUE))
    db.path.x <- unlist(strsplit(db.path, '/'))[length(unlist(strsplit(db.path, '/')))]
    db.path.x <- file.path("temp_vsearch", db.path.x)
    system2(vsearch.path, paste(" --usearch_global ", query.path, " --db ", db.path.x, " --blast6out ", b6.out, " --uc ", uc.out, " --iddef ",iddef, " --id 0.7 -query_cov 0.95 --maxaccepts 0 --maxrejects 0 --top_hits_only --strand both ", sep=""), stdout="", stderr="")
    #Clean up
    unlink(file.path(path, db.path.x),recursive=TRUE)
  }

  #:::::::::::::::::
  #Organize results
  #:::::::::::::::::
  message(cat(paste0("\033[97;", 40, "m","Determining closest species match", "\033[0m")))
  uc.results <- read.csv(uc.out, sep="\t", header = FALSE) %>%
    mutate(V4 = gsub("*", "51", .$V4, fixed=TRUE)) %>% mutate_at(vars(V4), as.numeric) %>%   # Fix percent ID column to make numeric
    mutate(V10 = gsub("*", "No_match", .$V10, fixed=TRUE)) %>%                               # Fix unmatched taxa column
    arrange(desc(V4)) %>% distinct(V9, .keep_all=TRUE) #%>% select(V1, V4, V9, V10) #%>% filter(V10 != "*")

  query.seqs <- Biostrings::readBStringSet(query.path)
  query.seqs.df <- as.data.frame(names(query.seqs)) %>%
    mutate(V9 = names(query.seqs), query_seq = paste(query.seqs)) %>%
    mutate(length= Biostrings::nchar(query_seq)) %>%
    mutate(Ns = stringr::str_count(.$query_seq, "[Nn]")) %>%
    select(-1)

  uc.results <- merge(uc.results, query.seqs.df, by="V9", all=TRUE)

  if(keep_temp_files==FALSE){
    #Clean up files
    unlink(file.path(path, uc.out),recursive=TRUE) #Remove UC file
    unlink(file.path(path, b6.out),recursive=TRUE) #Remove B6 file
    unlink(path,recursive=TRUE) #Remove temp directory
  }
  return(uc.results)
}
