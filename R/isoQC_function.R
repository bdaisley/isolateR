#' @name isoQC
#' @title Perform automated quality trimming of input .ab1 files
#' @description This function loads in ABIF files (.ab1 extension) and performs automatic quality trimming in batch mode.
#' @rdname isoQC
#' @export
#' @param input Path to directory with .ab1 files.
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
#' @seealso \code{\link{isoTAX}}, \code{\link{isoLIB}}
#' @return Returns quality trimmed Sanger sequences in FASTA format.
#' @importFrom stringr str_count
#' @importFrom seqinr write.fasta
#' @importFrom sangerseqR read.abif
#' @importFrom sangeranalyseR SangerRead
#' @examples
#' #Set path to directory containing example .ab1 files
#' fpath1 <- system.file("extdata/abif_examples/rocket_salad", package = "isolateR")
#'
#' #Step 1: Run isoQC function with default settings
#' isoQC.S4 <- isoQC(input=fpath1)
#'
#' #Show summary statistics
#' isoQC.S4

isoQC <- function(input=NULL,
                  export_html=TRUE,
                  export_csv=TRUE,
                  export_fasta=TRUE,
                  export_fasta_revcomp=FALSE,
                  verbose=FALSE,
                  exclude=NULL,
                  min_phred_score=20,
                  min_length=200,
                  sliding_window_cutoff=NULL,
                  sliding_window_size = 15,
                  date=NULL,
                  files_manual=NULL){

  
  # Reading files--------------------------------------------------------------------
  # Check input path
  if(is.null(input)) stop('Name of folder not supplied. Please check the execution script preamble to make sure it is entered correctly', call.=FALSE)
  if(!is.null(sliding_window_size)){
    if(sliding_window_size > 40){
      warning('"sliding_window_size" is set to above 40. It is recommended to set this parameter between 1-40 for best trimming results.')
    }
  }
  if(is.null(sliding_window_size)) {
    sliding_window_size = 15
    message(cat(paste0("\n", "\033[97;", 40, "m","Setting defaults for NULL inputs:", "\033[0m", "\033[35;", 49, "m"," 'sliding_window_size = 15'", "\033[0m" )))
  }
  
  #Setting folder paths for ABIF files----------------------------------------------------
  
  folder <- input
  setwd(folder)
  folder <- getwd()
  path <- stringr::str_replace_all(folder, '\\\\', '/')
  fname <- list.files(path)
  pattern <- 'ab1'
  abif_files <- stringr::str_subset(fname, pattern)
  if(is.null(files_manual)==FALSE){abif_files <- files_manual}
  
  #excluding specified files
  if(!is.null(exclude)) {
    
    # first check if files specified in exclude are in abif_files
    if(any(!exclude %in% abif_files)) {
      absent <- exclude[!exclude %in% abif_files]
      
      msg <- sprintf('%s could not be excluded because they were not found in your input folder. Only files in the input folder can be excluded.',
                     str_c(absent, collapse=', '))
      stop(msg, call.=FALSE)
    }
    else{
      abif_files <- abif_files[!abif_files %in% exclude]
      msg <- sprintf('%s have been excluded from your fasta file.',
                     str_c(exclude, collapse=', '))
      message(msg)
    }
  }
  
  #Checking overwrite inputs-------------------------------------------------
  if(!is.null(date)){
    date.formatted <- paste(str_pad(str_split_fixed(date, "_", 3)[,1], 2, pad = "0"),
                            str_pad(str_split_fixed(date, "_", 3)[,2], 2, pad = "0"),
                            str_pad(str_split_fixed(date, "_", 3)[,3], 2, pad = "0"), sep="_")
    message(cat(paste0("\n", "\033[97;", 40, "m","Setting date (YY/MM/DD) as: ",date.formatted, "\033[0m", "\n")))
  }
  
  if(!is.null(min_length)){
    min_length <- 200
  }
  
  #initializing objects--------------------------------------------------------------
  trim.start <- c()
  trim.end <- c()
  rawsparklist <- list()
  trimsparklist <- list()
  rawlist <- list()
  trimlist <- list()
  checkseq <- c()
  date.list <- c()
  length.list <- c()
  
  #set progress bar----------------------------------------------------------------
  message(cat(paste0("\033[97;", 95, "m","...Trimming ", paste(length(abif_files)) ," input files...", "\033[0m")))
  prog.bar.x <- txtProgressBar(min = 0,      			       # Minimum value of the progress bar
                               max = length(abif_files), # Maximum value of the progress bar
                               style = 3,    			       # Progress bar style (also available style = 1 and style = 2)
                               #width = 50,   			     # Progress bar width. Defaults to getOption("width")
                               char = "=")
  
  
  
  #Get max sequence length info-----------------------------------------------------------
  for(i in 1:length(abif_files)){
    fpath <- file.path(path, abif_files[i])
    length.list <- c(length.list, as.numeric(nchar(unlist((sangerseqR::read.abif(fpath))@data['PBAS.2'])))+5)
  }
  length.list.max <- max(length.list)
  
  for(i in 1:length(abif_files)){
    
    fpath <- file.path(path, abif_files[i])
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    # Basecalling steps
    #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
    abif.pre <- sangerseqR::read.abif(fpath)
    
    if(is.null(sliding_window_cutoff)){
      sliding_window_cutoff.x <- max(unlist(abif.pre@data['PCON.1']))*0.66
    } else {
      sliding_window_cutoff.x <- sliding_window_cutoff
    }
    
    if(!is.null(date)){
      date.list <- c(date.list, date.formatted)
    } else {
      date.list <- c(date.list, get_sanger_date(file=abif.pre))
    }
    
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Trimming file: ", fpath, "\033[0m")))}
    
    if(verbose==FALSE){
      invisible(capture.output(capture.output(abif1 <- tryCatch(
        {
          sangeranalyseR::SangerRead(readFeature           = "Forward Read",
                                     readFileName          = fpath,
                                     geneticCode           = GENETIC_CODE,
                                     TrimmingMethod        = "M2",
                                     M1TrimmingCutoff      = NULL,
                                     M2CutoffQualityScore  = round(sliding_window_cutoff.x),
                                     M2SlidingWindowSize   = sliding_window_size,
                                     baseNumPerRow         = 100,
                                     heightPerRow          = 200,
                                     signalRatioCutoff     = 0.33)
        },
        error = function(e) {
          message(cat(paste0("\033[97;", 41, "m","First trimming attempt failed, retrying with 'sliding window = 1'", "\033[0m")))
          sangeranalyseR::SangerRead(readFeature           = "Forward Read",
                                     readFileName          = fpath,
                                     geneticCode           = GENETIC_CODE,
                                     TrimmingMethod        = "M2",
                                     M1TrimmingCutoff      = NULL,
                                     M2CutoffQualityScore  = round(sliding_window_cutoff.x),
                                     M2SlidingWindowSize   = 1,
                                     baseNumPerRow         = 100,
                                     heightPerRow          = 200,
                                     signalRatioCutoff     = 0.33)
        }
      ), type="message"), type="output"))
    } else {abif1 <- tryCatch(
      {
        sangeranalyseR::SangerRead(readFeature           = "Forward Read",
                                   readFileName          = fpath,
                                   geneticCode           = GENETIC_CODE,
                                   TrimmingMethod        = "M2",
                                   M1TrimmingCutoff      = NULL,
                                   M2CutoffQualityScore  = round(sliding_window_cutoff.x),
                                   M2SlidingWindowSize   = sliding_window_size,
                                   baseNumPerRow         = 100,
                                   heightPerRow          = 200,
                                   signalRatioCutoff     = 0.33)
      },
      error = function(e) {
        message(cat(paste0("\033[97;", 41, "m","First trimming attempt failed, retrying with 'sliding window = 1'", "\033[0m")))
        sangeranalyseR::SangerRead(readFeature           = "Forward Read",
                                   readFileName          = fpath,
                                   geneticCode           = GENETIC_CODE,
                                   TrimmingMethod        = "M2",
                                   M1TrimmingCutoff      = NULL,
                                   M2CutoffQualityScore  = round(sliding_window_cutoff.x),
                                   M2SlidingWindowSize   = 1,
                                   baseNumPerRow         = 100,
                                   heightPerRow          = 200,
                                   signalRatioCutoff     = 0.33)
      }
    )}
    
    #start=abif1@QualityReport@trimmedStartPos
    abif1.start = NULL
    abif1.start <- abif1@QualityReport@trimmedStartPos + 5
    if(abif1.start >= abif1@QualityReport@rawSeqLength){abif1.start <- abif1@QualityReport@trimmedStartPos}
    if(abif1.start >= abif1@QualityReport@trimmedFinishPos){abif1.start <- abif1@QualityReport@trimmedStartPos}
    
    abif1.end <- nchar(paste(abif1@primarySeq))
    abif2 <- DNAStringSet(abif1@primarySeq, start=abif1.start, end=abif1.end)
    
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Searching for reverse primer endpoint: ", fpath, "\033[0m")))}
    
    if(verbose==FALSE){
      invisible(capture.output(capture.output(abif3 <- DECIPHER::TrimDNA(DNAStringSet(paste(abif2)),
                                                                         leftPatterns="",
                                                                         rightPatterns="AAAAAAAAAAAAAAAAAAAAAAAAA", #CTGCTGCCTYCCGTA
                                                                         minWidth=1,
                                                                         maxDistance = 0.2,
                                                                         minOverlap = 15,
                                                                         allowInternal=TRUE,
                                                                         threshold = 1,
                                                                         maxAverageError = 1,
                                                                         maxAmbiguities = 0.5,
                                                                         type="both"), type="message"), type="output"))
    } else {
      abif3 <- DECIPHER::TrimDNA(DNAStringSet(paste(abif2)),
                                 leftPatterns="",
                                 rightPatterns="AAAAAAAAAAAAAAAAAAAAAAAAA",
                                 minWidth=1,
                                 maxDistance = 0.2,
                                 minOverlap = 15,
                                 allowInternal=TRUE,
                                 threshold = 1,
                                 maxAverageError = 1,
                                 maxAmbiguities = 0.5,
                                 type="both")
    }
    
    #This code chunk checks if trimmed sequence region looks okay before proceeding to next step
    #----------------------------------------------------------------------------------------------
    #if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Removing anything beyond reverse primer: ", fpath, "\033[0m")))} #Not implemented yet, estimating 30mer primer length
    abif3.end <- end(abif3[[1]]) - 30  #+ nchar("CTGCTGCCTYCCGTA")
    if(abif3.end > width(abif3[[1]])){abif3.end <- width(abif3[[1]])}
    if(width(abif3[[1]]) <= 300){abif3.end <- nchar(paste(abif2)) }
    if(abif3.end > (abif1@QualityReport@trimmedFinishPos - abif1.start)){abif3.end <- (abif1@QualityReport@trimmedFinishPos - abif1.start)}
    if(abif3.end > length.list.max){abif3.end <- length.list.max}
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Trimming error check: 1 (Passed)", "\033[0m")))}
    abif1.start.new <- abif1.start
    if(abif1.start <= 50 & abif3.end >= 51){abif1.start.new <- 50}
    abif <- DNAStringSet(abif1@primarySeq, start=abif1.start.new, end = (abif1.start + abif3.end -1))
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Trimming error check: 2 (Passed)", "\033[0m")))}
    abif.scores.xx <- abif1@QualityReport@qualityPhredScores[abif1.start.new:(abif1.start + abif3.end -1)]
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Trimming error check: 3 (Passed)", "\033[0m")))}
    abif.scores <- c(rep(0.1, abif1.start.new), abif.scores.xx, rep(0.1, length.list.max-(((abif1.start + abif3.end) -1))))
    trim.start <- c(trim.start,abif1.start.new)
    trim.end <- c(trim.end,(abif1.start + abif3.end -1))
    rawseq <- as.character(paste(abif))
    
    if(verbose==TRUE){message(cat(paste0("\033[97;", 40, "m","Saving seq info: ", fpath, "\033[0m")))}
    
    if(verbose==TRUE){
      msg <- sprintf('Processing "%s"', abif_files[i])
      message(msg)
    }
    
    
    #:::::::::::::::::::::::::::::::
    # Save seq + quality scores
    #:::::::::::::::::::::::::::::::
    rawspark.adj <- abif1@QualityReport@qualityPhredScores
    if(length(abif1@QualityReport@qualityPhredScores) >= length.list.max){rawspark.adj <- rawspark.adj[1:length.list.max]}
    rawsparklist[[i]] <- c(rawspark.adj, rep(0.1, length.list.max-length(rawspark.adj)))
    trimsparklist[[i]] <- abif.scores
    rawlist[[i]] <- paste(abif1@primarySeqRaw)
    trimlist[[i]] <- paste(abif)
    
    # building trim check for export-------------------------------------------------
    entry <- data.frame(date=date.list[i],
                        filename=abif_files[i],
                        trim.start.pos = abif1.start.new,
                        trim.end.pos = abif1.start + abif3.end -1,
                        #spark_data = list(sparklist[[i]]),
                        phred_spark_raw = rbind(list(rawsparklist[[i]])),
                        phred_raw = round(mean(rawsparklist[[i]]), 2),
                        phred_trim = round(mean(abif.scores.xx),2 ),
                        length_raw = nchar(rawlist[[i]]),
                        length_trim = nchar(trimlist[[i]]),
                        Ns_raw = stringr::str_count(rawlist[[i]], "[Nn]"),
                        Ns_trim = stringr::str_count(trimlist[[i]], "[Nn]"),
                        seqs_raw = rawlist[[i]],
                        seqs_trim = trimlist[[i]]) %>%
      mutate(decision = ifelse(.$length_trim < min_length | .$phred_trim < min_phred_score, "Fail", "Pass"))
    
    # building trim check for export-------------------------------------------------
    checkseq <- rbind(checkseq, entry)
    if(verbose==FALSE){setTxtProgressBar(prog.bar.x, i)}
    Sys.sleep(0.01)
  }
  # end of file loop
  Sys.sleep(0.02)
  
  
  # building isoQC file--------------------------------------------------------
  
  isoQC <- new("isoQC", input=paste(path))
  
  isoQC@date <- checkseq$date
  isoQC@trim.start.pos <- checkseq$trim.start.pos
  isoQC@trim.end.pos <- checkseq$trim.end.pos
  isoQC@filename <- checkseq$filename
  isoQC@phred_spark_raw <- checkseq$phred_spark_raw
  isoQC@phred_raw <- checkseq$phred_raw
  isoQC@phred_trim <- checkseq$phred_trim
  isoQC@length_raw <- checkseq$length_raw
  isoQC@length_trim <- checkseq$length_trim
  isoQC@Ns_raw <- checkseq$Ns_raw
  isoQC@Ns_trim <- checkseq$Ns_trim
  isoQC@seqs_raw <- checkseq$seqs_raw
  isoQC@seqs_trim <- checkseq$seqs_trim
  isoQC@decision <- checkseq$decision
  
  #:::::::::::::::::::::::::::::::
  # Output results
  #:::::::::::::::::::::::::::::::
  isoQC.df <- S4_to_dataframe(isoQC) %>%
    select(date,
           filename,
           seqs_raw,
           phred_raw,
           Ns_raw,
           length_raw,
           phred_spark_raw,
           #arrcol,
           seqs_trim,
           phred_trim,
           Ns_trim,
           length_trim,
           decision)
  
  #building output file--------------------------------------------------------------
  
  seq.warnings <- (isoQC.df %>% filter(decision!="Pass"))$filename
  
  if(is.null(seq.warnings)==FALSE){
    decision.txt <- paste(seq.warnings, collapse="\r\n     ")
    message(cat(paste0("\033[97;", 40, "m","\r\nThe following sequences failed with a length of <", min_length," bp and/or a quality score <", min_phred_score ," after trimming:","\033[0m","\n     ",
                       "\033[0;", 95, "m", decision.txt,"\033[0m","\n")))
  }
  
  message(cat(paste0("\n", "\033[97;", 40, "m","Export directory:", "\033[0m",
                     "\033[0;", 32, "m", " ", file.path(path, "isolateR_output"), "\033[0m","\n")))
  
  suppressWarnings(dir.create(file.path(path, "isolateR_output")))
  fname <- file.path(path,"isolateR_output", paste0("01_isoQC_results", sep=""))
  fname_pass <- file.path(path,"isolateR_output", paste0("01_isoQC_trimmed_sequences_PASS", sep=""))
  fname_fail <- file.path(path,"isolateR_output", paste0("01_isoQC_trimmed_sequences_FAIL", sep=""))
  isoQC.df.pass <- isoQC.df %>% filter(!filename %in% seq.warnings)
  isoQC.df.fail <- isoQC.df %>% filter(filename %in% seq.warnings)
  
  #export HTML file----------------------------------------------------------------------
  if(export_html == TRUE) {
    fname_html <- paste0(fname, ".html", sep="")
    export_html(isoQC)
    message(cat(paste0("\033[97;", 40, "m","HTML results exported: ", "\033[0m",
                       "\033[0;", 32, "m", " ", file.path(path, "isolateR_output", unlist(strsplit(fname_html, '/'))[length(unlist(strsplit(fname_html, '/')))]),"\033[0m", "\n")))
  }
  
  #export CSV files----------------------------------------------------------------------
  if(export_csv == TRUE) {
    #PASS sequences
    #---------------
    fname_csv_pass <- paste0(fname_pass, ".csv", sep="")
    if(isEmpty(isoQC.df.pass[1])){
      message(cat(paste0("\033[97;", 40, "m","CSV file with [PASS] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", "No sequences failed, file not exported.","\033[0m")))
    } else {
      write.csv(file=fname_csv_pass,isoQC.df.pass, row.names = FALSE)
      message(cat(paste0("\033[97;", 40, "m","CSV results exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, "isolateR_output", unlist(strsplit(fname_csv_pass, '/'))[length(unlist(strsplit(fname_csv_pass, '/')))]),"\033[0m",
                         "\033[0;", 31, "m", "  <--- Required in Step 2: 'isoTAX'","\033[0m")))
    }
    #FAIL sequences
    #---------------
    fname_csv_fail <- paste0(fname_fail, ".csv", sep="")
    if(isEmpty(isoQC.df.fail[1])){
      message(cat(paste0("\033[97;", 40, "m","CSV file with [FAIL] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", "No sequences failed, file not exported.","\033[0m")))
    } else {
      write.csv(file=fname_csv_fail,isoQC.df.fail, row.names = FALSE)
      message(cat(paste0("\033[97;", 40, "m","CSV results exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, "isolateR_output", unlist(strsplit(fname_csv_fail, '/'))[length(unlist(strsplit(fname_csv_fail, '/')))]),"\033[0m")))
    }
  }
  
  
  #export trimmed sequences in FASTA file-------------------------------------------------
  if(export_fasta == TRUE) {
    #PASS sequences
    #---------------
    fname_fasta_pass <- paste0(fname_pass, ".fasta", sep="")
    if(isEmpty(isoQC.df.pass[1])){
      message(cat(paste0("\033[97;", 40, "m","FASTA file with [PASS] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", "No sequences failed, file not exported.","\033[0m")))
    } else {
      trimlist_pass <- as.list((unlist(trimlist))[which(!isoQC.df$filename %in% seq.warnings)])
      abif_files_pass <- abif_files[!abif_files %in% seq.warnings]
      seqinr::write.fasta(sequences=trimlist_pass, names=abif_files_pass, as.string=TRUE, nbchar = 1000,
                  file.out=fname_fasta_pass)
      message(cat(paste0("\033[97;", 40, "m","FASTA file with [PASS] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, "isolateR_output", unlist(strsplit(fname_fasta_pass, '/'))[length(unlist(strsplit(fname_fasta_pass, '/')))]),"\033[0m")))
    }
    #FAIL sequences
    #---------------
    fname_fasta_fail <- paste0(fname_fail, ".fasta", sep="")
    if(isEmpty(isoQC.df.fail[1])){
      message(cat(paste0("\033[97;", 40, "m","FASTA file with [FAIL] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", "No sequences failed, file not exported.","\033[0m")))
    } else {
      trimlist_fail <- as.list((unlist(trimlist))[which(isoQC.df$filename %in% seq.warnings)])
      abif_files_fail <- abif_files[abif_files %in% seq.warnings]
      seqinr::write.fasta(sequences=trimlist_fail, names=abif_files_fail, as.string=TRUE, nbchar = 1000,
                  file.out=fname_fasta_fail)
      message(cat(paste0("\033[97;", 40, "m","FASTA file with [FAIL] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, "isolateR_output", unlist(strsplit(fname_fasta_fail, '/'))[length(unlist(strsplit(fname_fasta_fail, '/')))]),"\033[0m", "\n")))
    }
  }
  if(export_fasta_revcomp==TRUE) {
    #PASS sequences
    #---------------
    fname_fasta_revcomp_pass <- paste0(fname_pass, "_revcomp.fasta", sep="")
    if(isEmpty(isoQC.df.pass[1])){
      message(cat(paste0("\033[97;", 40, "m","FASTA (reverse complement) file with [PASS] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", "No sequences failed, file not exported.","\033[0m")))
    } else {
      trimlist_pass <- as.list(paste(Biostrings::reverseComplement(DNAStringSet((unlist(trimlist))[which(!isoQC.df$filename %in% seq.warnings)]))))
      abif_files_pass <- abif_files[!abif_files %in% seq.warnings]
      suppressWarnings({
        seqinr::write.fasta(sequences=trimlist_pass, names=abif_files_pass, as.string=TRUE, nbchar = 1000,
                    file.out=fname_fasta_revcomp_pass)
      })
      message(cat(paste0("\033[97;", 40, "m","FASTA (reverse complement) file with [PASS] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, "isolateR_output", unlist(strsplit(fname_fasta_revcomp_pass, '/'))[length(unlist(strsplit(fname_fasta_revcomp_pass, '/')))]),"\033[0m")))
    }
    #FAIL sequences
    #---------------
    fname_fasta_revcomp_fail <- paste0(fname_fail, "_revcomp.fasta", sep="")
    if(isEmpty(isoQC.df.fail[1])){
      message(cat(paste0("\033[97;", 40, "m","FASTA (reverse complement) file with [FAIL] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", "No sequences failed, file not exported.","\033[0m")))
    } else {
      trimlist_fail <- as.list(paste(Biostrings::reverseComplement(DNAStringSet((unlist(trimlist))[which(isoQC.df$filename %in% seq.warnings)]))))
      abif_files_fail <- abif_files[abif_files %in% seq.warnings]
      suppressWarnings({
        seqinr::write.fasta(sequences=trimlist_fail, names=abif_files_fail, as.string=TRUE, nbchar = 1000,
                    file.out=fname_fasta_revcomp_fail)
      })
      message(cat(paste0("\033[97;", 40, "m","FASTA (reverse complement) file with [FAIL] sequences exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, "isolateR_output", unlist(strsplit(fname_fasta_revcomp_fail, '/'))[length(unlist(strsplit(fname_fasta_revcomp_fail, '/')))]),"\033[0m")))
    }
  }
  
  #remove the folder of misc files for the script
  unlink(paste0(fname,"_files/"),recursive = TRUE)
  return(isoQC)
}

