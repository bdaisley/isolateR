#' @name sanger_consensus
#' @title Overlap multiple paired Sanger sequences in batch.
#' @description This function loads in the CSV results table from isoQC and merges related sequences based on user input. Original file 
#' names before isoQC step need to have a common prefix and differentiating suffixes. (e.g. SAMPLE_01_F.ab1, SAMPLE_01_R.ab1). After
#' aligning paired sequences, the consensus sequence is extracted and priority is given to the read with higher quality. Phred quality
#' scores are reassigned in the final output table in a basic way by taking the mean of both input sequences.
#'
#' Note: This function is designed to be used after the isoQC step and before the isoTAX step.
#' @rdname sanger_consensus
#' @export
#' @param input Path of CSV output file from isoQC step.
#' @param suffix Regex-friendly suffix for denoting filename groupings. Default="_F.ab1|_R.ab1" for the common scenario of Sanger sequencing a marker gene in forward and reverse. 
#' Direction of sequences including reverse complements will be automatically detected.
#' @seealso \code{\link{isoQC}}, \code{\link{isoTAX}}
#' @return Returns merged pairs of Sanger sequences in FASTA format.
#' @importFrom DECIPHER AlignSeqs
#' @importFrom DECIPHER ConsensusSequence
#' @importFrom Biostrings reverseComplement
#' @importFrom stringr str_subset
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_split_fixed
#' @importFrom stringr str_remove
#' 


sanger_consensus <- function(input=NULL,
                             suffix="_F.ab1|_R.ab1"){
  
  #::::::::::::::::::
  #---Reading files
  #::::::::::::::::::
  # Check input path
  if(is.null(input)) stop('Name of folder not supplied. Please check the execution script preamble to make sure it is entered correctly', call.=FALSE)
  
  #::::::::::::::::::
  #---Input files
  #::::::::::::::::::
  input.df <- read.csv(input) %>% mutate(no_suffix= gsub(suffix, "",filename)) %>% mutate(row_counter=seq_along(1:n()))
  
  paired.seqs <- input.df %>% dplyr::group_by(no_suffix) %>% dplyr::filter(!n() < 2) %>% dplyr::ungroup()
  unpaired.seqs <- input.df %>% dplyr::group_by(no_suffix) %>% dplyr::filter(n() < 2) %>% dplyr::ungroup()
  
  if(length(unique(paired.seqs$no_suffix)) == 0 ) stop("There were 0 groups found based on common part of filename. Try adjusting 'suffix' parameter.")
  
  message("Detected ", length(unique(paired.seqs$no_suffix)), " unique group(s) with suffix provided.")
  
  if(length(unique(unpaired.seqs$no_suffix)) != 0){message("Warning: ", length(unique(unpaired.seqs$no_suffix)), " files had no matching pair after trimming suffix.",
                                                           " The file(s) will be skipped during alignment steps but will still be included in final output.")}

  #::::::::::::::::::
  #---Align seqs
  #::::::::::::::::::
  seq.groupings <- unique(paired.seqs$no_suffix)
  
  collector.list <- lapply(1:length(seq.groupings), function(i){
    #Filter df to get only group of interest
    paired.seqs.filt <- paired.seqs %>% filter(no_suffix == seq.groupings[i]) # ************change back to i
    
    #---------------Set initial values
    seqs <- Biostrings::DNAStringSet(paired.seqs.filt$seqs_trim)
    seqs.rc <- reverseComplement(seqs)
    seqs.qual <- lapply(1:nrow(paired.seqs.filt), function(x) rep(paired.seqs.filt$phred_trim[x], paired.seqs.filt$length_trim[x]))
    seqs.qual.rc <- rev(seqs.qual)
    
    init.seq <- seqs[1]
    init.qual <- seqs.qual[[1]]
    failed.aln <- c()
    #---------------
    while(length(seqs) > 0){
      #Evaluate input sequences as is
      eval.1 <-unlist(lapply(seqs, function(x){Biostrings::pairwiseAlignment(pattern=x, subject=init.seq, type = "overlap", substitutionMatrix = NULL, scoreOnly=TRUE)}))
      #Evaluate reverse complement of input sequences
      eval.2 <-unlist(lapply(seqs.rc, function(x){Biostrings::pairwiseAlignment(pattern=x, subject=init.seq, type = "overlap", substitutionMatrix = NULL, scoreOnly=TRUE)}))
      #Get best matching seq based on alignment score
      best.match.index <- ifelse(max(eval.1) > max(eval.2), which.max(eval.1), which.max(eval.2))
      best.match.score <- ifelse(max(eval.1) > max(eval.2), eval.1[which.max(eval.1)], eval.2[which.max(eval.2)])
      best.match <- ifelse(max(eval.1) > max(eval.2), seqs[which.max(eval.1)], seqs.rc[which.max(eval.2)])
      best.match.qual <- ifelse(max(eval.1) > max(eval.2), seqs.qual[[which.max(eval.1)]], seqs.qual.rc[[which.max(eval.2)]])
      
      #DECIPHER functions------------
      eval.3 <- DECIPHER::AlignSeqs(DNAStringSet(c(init.seq, best.match, gsub("-", "", alignedPattern(Biostrings::pairwiseAlignment(pattern=best.match,subject=init.seq,type = "overlap", substitutionMatrix = NULL))))), verbose=FALSE)
      
      #---Debugging: Check alignment
      if(best.match.score < 20) message("***Alignment failed for ", unique(paired.seqs.filt$no_suffix),": Alignment score (", round(best.match.score, 2) ,
                                        ") less than 20 for at least one pair of sequences. This could be due to very short overlaps, or poor sequence quality. ",
                                        "Manual inspection of isoQC output is recommended. Try changing isoQC parameters 'sliding_window_cutoff' and 'sliding_window_size'")

      #DECIPHER::BrowseSeqs(eval.3)
      #-----------------------------
      #Get consensus seq
      seqs.aln.con <- DECIPHER::ConsensusSequence(eval.3, threshold = 0.05, ambiguity = TRUE)
      #Combine query seqs + consensus side-by-side in dataframe
      seqs.aln.con.manual <- as.data.frame(rbind(c(unlist(as.character(eval.3)),
                                                   as.character(seqs.aln.con))))
      #Split sequences into individual nucleotides
      seqs.aln.con.manual.l <- lapply(1:ncol(seqs.aln.con.manual), function(x){
        unlist.out1 <- unlist(strsplit(seqs.aln.con.manual[,x], "(?<=[-A-Z])", perl=TRUE))
        unlist.out <- unlist(strsplit(unlist.out1, "(?<=[\\+])", perl=TRUE))
        return(unlist.out)})
      
      #Bind rows of list
      seqs.aln.con.manual <- suppressMessages(bind_cols(seqs.aln.con.manual.l)) %>% `colnames<-`(c("seq1","seq2", "seq3", "auto")) %>% 
        mutate(qual_1=rep(unique(init.qual))) %>% mutate(qual_2=rep(unique(best.match.qual))) %>%
        mutate(consensus_seq = ifelse((seq1 == "-" | seq2 == "-"), auto, ifelse(qual_1 > qual_2, seq1, seq2))) %>%
        mutate(consensus_seq = ifelse(consensus_seq=="+", ifelse(qual_1 > qual_2, seq1, seq2), consensus_seq)) %>%
        filter(consensus_seq!="-")
      
      #Update index references
      init.seq <- paste(seqs.aln.con.manual$consensus_seq, collapse="")
      init.qual <- round(mean(c(seqs.aln.con.manual$qual_1, seqs.aln.con.manual$qual_2)), 2)
      
      seqs <- seqs[-best.match.index]
      seqs.rc <- seqs.rc[-best.match.index]
      seqs.qual <- seqs.qual[-best.match.index]
      seqs.qual.rc <- seqs.qual.rc[-best.match.index]
      #seqs.qual.rc2 <- sapply(seqs.qual.rc, base::toString)
      
      if(best.match.score <20){
        seqs <- NULL
        failed.aln <- c(failed.aln, unique(paired.seqs.filt$no_suffix))
        }
      
    } #end of while loop
    
    if(length(failed.aln) > 0){
      seqs.aln.con.df <- paired.seqs.filt %>% 
        mutate(phred_trim = round(mean(phred_trim), 2)) %>% 
        mutate(decision="Pass") %>%
        .[1,] %>%
        mutate(seqs_trim = init.seq) %>%
        mutate(Ns_trim=nchar(as.character(seqs_trim)) - stringr::str_count(as.character(seqs_trim), "[ACTG]")) %>%
        mutate(length_trim=nchar(as.character(seqs_trim))) %>%
        mutate(aln=paste("fail"))
    } else {
      seqs.aln.con.df <- paired.seqs.filt %>% 
        mutate(phred_trim = round(mean(phred_trim), 2)) %>% 
        mutate(decision="Pass") %>%
        .[1,] %>%
        mutate(seqs_trim = init.seq) %>%
        mutate(Ns_trim=nchar(as.character(seqs_trim)) - stringr::str_count(as.character(seqs_trim), "[ACTG]")) %>%
        mutate(length_trim=nchar(as.character(seqs_trim))) %>%
        mutate(aln=paste("pass"))
    }
    return(seqs.aln.con.df) 
  }) #end of collector.list function

  #Bind all overlapped sequences back together and add back in unpaired seqs
  if(any((bind_rows(collector.list))$aln=="pass")){
    if(length(paste(unique((bind_rows(collector.list) %>% filter(aln=="fail"))$no_suffix), sep="|"))==0){
      collector.list.df <- bind_rows(collector.list) %>%
        {tmp <<- .} %>%
        filter(aln=="pass") %>%
        mutate(seqs_raw=rep("xNA")) %>%
        mutate(phred_raw=rep(0)) %>%
        mutate(Ns_raw=rep(0)) %>%
        mutate(length_raw=rep(0)) %>%
        mutate(phred_spark_raw=rep(c("0_0_0"))) %>%
        select(-aln) %>%
        mutate(filename= paste(no_suffix, "_consensus.",stringr::str_split_fixed(.$filename[1], "[.]", 2)[,2], sep="")) %>%
        rbind(., 
              unpaired.seqs) %>%
        select(-no_suffix, -row_counter)
    } else {
      collector.list.df <- bind_rows(collector.list) %>%
        {tmp <<- .} %>%
        filter(aln=="pass") %>%
        mutate(seqs_raw=rep("xNA")) %>%
        mutate(phred_raw=rep(0)) %>%
        mutate(Ns_raw=rep(0)) %>%
        mutate(length_raw=rep(0)) %>%
        mutate(phred_spark_raw=rep("NA")) %>%
        select(-aln) %>%
        mutate(filename= paste(no_suffix, "_consensus.",stringr::str_split_fixed(.$filename[1], "[.]", 2)[,2], sep="")) %>%
        rbind(., 
              (paired.seqs %>% filter(grepl(paste(unique((tmp %>% filter(aln=="fail"))$no_suffix), sep="|"), no_suffix))),
              unpaired.seqs) %>%
        select(-no_suffix, -row_counter)
    }
  } else{
    collector.list.df <- rbind(paired.seqs, 
                               unpaired.seqs) %>% 
      select(-no_suffix, -row_counter)
  }
  

  #:::::::::::::::::::::::
  #---Export CSV files
  #:::::::::::::::::::::::
  csv.out.path <- paste(stringr::str_split_fixed(input, "[.]csv", 2)[,1], "_consensus.csv", sep="")
  write.csv(collector.list.df, csv.out.path, row.names = FALSE)
  
  message(cat(paste0("\033[97;", 40, "m","CSV file with merged sequences exported:", "\033[0m",
                     "\033[0;", 32, "m", " ", csv.out.path,"\033[0m",
                     "\033[0;", 31, "m", "  <--- Input for Step 2: 'isoTAX'","\033[0m")))
  
  #:::::::::::::::::::::::::::::::::::::::::::
  #---Return data as isoQC class object in R
  #:::::::::::::::::::::::::::::::::::::::::::
  isoQC <- new("isoQC", input=paste(getwd()))
  isoQC@date <- collector.list.df$date
  isoQC@trim.start.pos <- rep(0, nrow(collector.list.df))
  isoQC@trim.end.pos <- rep(0, nrow(collector.list.df))
  isoQC@filename <- collector.list.df$filename
  isoQC@phred_spark_raw <- as.list(collector.list.df$phred_spark_raw)
  isoQC@phred_raw <- collector.list.df$phred_raw
  isoQC@phred_trim <- collector.list.df$phred_trim
  isoQC@length_raw <- collector.list.df$length_raw
  isoQC@length_trim <- collector.list.df$length_trim
  isoQC@Ns_raw <- collector.list.df$Ns_raw
  isoQC@Ns_trim <- collector.list.df$Ns_trim
  isoQC@seqs_raw <- collector.list.df$seqs_raw
  isoQC@seqs_trim <- collector.list.df$seqs_trim
  isoQC@decision <- collector.list.df$decision
  
  return(isoQC)
}
