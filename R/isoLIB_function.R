#' @title Generate new strain library or add to existing one.
#' @name isoLIB
#' @rdname isoLIB
#' @description This function creates a strain library by grouping closely related strains of interest based on sequence similarity.
#' For adding new sequences to an already-established strain library, specify the .CSV file path of the older strain library using the 'old_lib_csv" parameter.
#' @export
#' @param input Path of CSV output file from isoTAX step.
#' @param old_lib_csv Optional: Path of CSV output isoLIB file or combined isoLIB file from previous run(s)
#' @param method Method used for grouping sequences. Either 1) "dark_mode", or 2) "closest_species" (Default="dark_mode").
#' - Method 1 ("dark_mode") performs agglomerative hierarchical-based clustering to group similar sequences based on pairwise identity (see 'id' parameter)
#' and then within each group, attempts to assign the longest sequence with the most top hits as the group representative. This method is tailored for capturing distinct 
#' strains which may represent novel taxa (i.e. microbial dark matter) during isolation workflows. As such, the sequence representatives chosen in each group will not always 
#' have the highest % identity to the closest matching type strain. In some cases, sequence members within a group may also have different taxonomic classifications
#' due to them having close to equidistant % identities to different matching type strain material – indicative of a potentially novel taxonomic grouping.
#' - Method 2 ("closest_species") groups similar sequences based on their closest matching type strain. For each unique grouping, this results in all sequence members having
#' the same taxonomic classification. The longest sequence with the highest % identity to the closest matching type strain will be assigned as the group representative.
#' Note: The "id" parameter is only used for Method 1 ("dark_mode") and otherwise ignored if using Method 2 ("closest_species").
#' @param group_cutoff (Default=0.995) Similarity threshold based on pairwise identity (0-1) for delineating between sequence groups. 1 = 100% identical/0.995=0.5% difference/0.95=5.0% difference/etc.
#' Used only if method="dark_mode", otherwise ignored.
#' @param keep_old_reps (Default=TRUE) If TRUE, original sequence representatives from old library will be kept when merging with new library. 
#' If FALSE, sequence group representatives will be recalculated after combining old and new libraries.
#' Ignored if old_lib_csv=NULL.
#' @param export_html (Default=TRUE) Output the results as an HTML file
#' @param export_csv (Default=TRUE) Output the results as a CSV file.
#' @param include_warnings (Default=FALSE) Whether or not to keep sequences with poor alignment warnings from Step 2 'isoTAX' function. Set TRUE to keep warning sequences, and FALSE to remove warning sequences.
#' @param vsearch_path Path of VSEARCH software if manually downloaded in a custom directory. If NULL (Default), will attempt automatic download.
#' @param phylum_threshold Percent sequence similarity threshold for phylum rank demarcation
#' @param class_threshold Percent sequence similarity threshold for class rank demarcation
#' @param order_threshold Percent sequence similarity threshold for order rank demarcation
#' @param family_threshold Percent sequence similarity threshold for family rank demarcation
#' @param genus_threshold Percent sequence similarity threshold for genus rank demarcation
#' @param species_threshold Percent sequence similarity threshold for species rank demarcation
#' @seealso \code{\link{isoTAX}}, \code{\link{isoLIB}}
#' @return Returns an isoLIB class object. Default taxonomic cutoffs for phylum (75.0), class (78.5), order (82.0), family (86.5), genus (94.5), and species (98.7) demarcation are based on Yarza et al. 2014, Nature Reviews Microbiology (DOI:10.1038/nrmicro3330)
#' @importFrom utils download.file
#' @importFrom utils untar
#' @importFrom R.utils gunzip
#' @importFrom stringr str_sub
#' @examples
#' #Set path to directory containing example .ab1 files
#' fpath1 <- system.file("extdata/abif_examples/rocket_salad", package = "isolateR")
#'
#' #Step 1: Run isoQC function with default settings
#' isoQC.S4 <- isoQC(input=fpath1)
#'
#' #Step 2: Run isoTAX function with default settings
#' fpath2 <- file.path(fpath1, "isolateR_output/01_isoQC_trimmed_sequences_PASS.csv")
#' isoTAX.S4 <- isoTAX(input=fpath2)
#'
#' #Step 3: Run isoLIB function with default settings
#' fpath3 <- file.path(fpath1, "isolateR_output/02_isoTAX_results.csv")
#' isoLIB.S4 <- isoLIB(input=fpath3)
#'
#' #Show summary statistics
#' isoLIB.S4

isoLIB <- function(input=NULL,
                   old_lib_csv=NULL,
                   method="dark_mode",
                   group_cutoff = 0.995,
                   keep_old_reps=TRUE,
                   export_html=TRUE,
                   export_csv=TRUE,
                   include_warnings=TRUE,
                   vsearch_path=NULL,
                   phylum_threshold=75.0,
                   class_threshold=78.5,
                   order_threshold=82.0,
                   family_threshold=86.5,
                   genus_threshold=94.5,
                   species_threshold=98.7){
  
  #Input checks------------------------------------------------------------
  #Set alternatives
  method <- gsub("^ml$|^ML$", "maximum_likelihood", method)
  method <- gsub("^cs$|^CS$|^closest$|^species$|^s$", "closest_species", method)
  method <- gsub("^dm$|^DM$|^dark$|^darkmode$|^d$", "dark_mode", method)
  if(!grepl("^maximum_likelihood$|^closest_species$|^dark_mode$", method)) stop('Mode selection incorrect. Please choose one of: "maximum_likelihood", "closest_species", or "dark_mode"')
  if(!is.null(group_cutoff)){
    if(!is.numeric(group_cutoff)) stop("'group_cutoff' is not numeric. Set between 0-1 (1=no difference i.e., identical sequences)", call.=FALSE)
    if(group_cutoff <0 | group_cutoff >1) stop("Wrong 'group_cutoff' format. Set between 0-1 (1=no difference i.e., identical sequences)", call.=FALSE)
  }
  
  #Set paths------------------------------------------------------------
  
  new_lib <- stringr::str_replace_all(input, '\\\\', '/')
  new_lib_path <- paste(unlist(strsplit(new_lib, '/'))[1:(length(unlist(strsplit(new_lib, '/')))-1)],collapse="/")
  new_lib_file <- read.csv(new_lib) #, row.names = 1) removing row.names setting
  if(nrow(new_lib_file) < 2) stop("Command terminated. isoLIB requires at least 2 sequences to generate a library.", call.=FALSE)
  
  setwd(new_lib_path)
  
  #-------------------------------------------------
  #:::::::::::::::::::::::::::
  #Download VSEARCH software
  #:::::::::::::::::::::::::::
  
  path = new_lib_path
  vsearch.path.dl <- file.path(system.file(package="isolateR"), "vsearch")
  suppressWarnings(dir.create(vsearch.path.dl))
  vsearch_files <- stringr::str_subset(dir(vsearch.path.dl, full.names = FALSE), 'vsearch')
  
  if(is.null(vsearch_path)){
    vsearch.path <- get_vsearch()
  } else {
    vsearch.path <- vsearch_path
  }

  #::::::::::::
  #Stage paths
  #::::::::::::
  message(cat(paste0("\n", "\033[0;", 32, "m","Staging files for de-replication", "\033[0m", "\n")))
  
  temp_csv <- "03_isoLIB_temp.csv"
  temp_fasta <- "03_isoLIB_temp.fasta"
  temp_drep <- paste0("03_isoLIB_temp.drep", sep="")
  
  
  #Make unique strain names incase any duplicated-----------------------
  
  new_lib_file <- new_lib_file %>% arrange(desc(length_trim))
  
  #Check for duplicate names
  new_lib_file <- new_lib_file %>% 
    mutate(undup = stringr::str_split_fixed(filename, "_dup", 2)[,1]) %>% 
    mutate(duplicated = duplicated(undup)) %>%
    mutate(duplicated2 = paste(undup, duplicated, sep="___")) %>%
    group_by(duplicated2) %>%
    mutate(filename = ifelse(any(duplicated), paste(undup, "_dup", rep(1:n()), sep=""), filename)) %>%
    ungroup()
  
  if(length((new_lib_file %>% filter(duplicated==TRUE))$filename) > 0){
    message(cat(paste0("\033[0;", 31, "m", "Duplicated filenames detected: ", paste(unique((new_lib_file %>% filter(duplicated==TRUE))$undup), collapse="|"), "\033[0m")))
    message(cat(paste0("\033[0;", 30, "m", "Renaming duplicates with suffix '_dup1','_dup2', '_dup3', etc.", "\033[0m")))
    message(cat(paste0("\033[0;", 30, "m",  paste((new_lib_file %>% filter(duplicated==TRUE))$filename, collapse="\n"), "\033[0m")))
  }
  
  new_lib_file <- new_lib_file %>% select(-undup, -duplicated, -duplicated2)
  
  new_lib_file$unique_filename <- paste(stringr::str_pad(1:nrow(new_lib_file), 4, pad="0"), new_lib_file$filename, sep="_")
  
  
  ###########################################################################################################
  ###########################################################################################################
  ########################################################################################################### 3 'method' options
  ###########################################################################################################
  ###########################################################################################################
  
  ############################################################################
  ############################################################################
  ###############                                              ###############
  ###############  START /// method = "maximum_likelihood" /// ###############
  ###############                                              ###############
  ############################################################################
  ############################################################################
  
  if(method=="maximum_likelihood"){
    new.seqs <- setNames(new_lib_file$seqs_trim, new_lib_file$unique_filename)
    #Align sequences
    cluster.aln <- DECIPHER::AlignSeqs(DNAStringSet(new.seqs))
    #Calculate distance matrix
    cluster.dist <- DECIPHER::DistanceMatrix(cluster.aln, type="matrix")
    #Cluster sequences
    cluster.list <- DECIPHER::TreeLine(myXStringSet = cluster.aln,
                                       myDistMatrix = cluster.dist,
                                       cutoff = 1-group_cutoff,
                                       method = "ML", #UPGMA
                                       model="JC69",
                                       maxTime=0.02,
                                       processors=NULL,
                                       type="both")
    
    merged.drep1 <- new_lib_file %>% 
      mutate(sequence_group.x = cluster.list[[1]]$cluster) %>%
      mutate(row_no = row_number()) %>% 
      arrange(Ns_trim) %>%
      arrange(desc(phred_trim)) %>%
      arrange(desc(length_trim)) %>%
      group_by(sequence_group.x) %>%
      mutate(sequence_group = first(unique_filename)) %>%
      ungroup() %>%
      mutate(sequence_group = sapply(1:length(.$sequence_group),function(x) stringr::str_sub(.$sequence_group[x], 6, nchar(.$sequence_group[x])))) %>%
      mutate(representative = ifelse(sequence_group == filename, "yes", "no")) %>%
      arrange(row_no) %>%
      dplyr::select(sequence_group, date, filename, representative, seqs_trim, phred_trim, Ns_trim, length_trim, closest_match, NCBI_acc, ID,
                    rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species)
    
    #---Sanity check
    message(cat(paste0("\033[0;", 32, "m","A total of ", length(unique(merged.drep1$sequence_group)),
                       " unique sequence groups detected based on natural clustering on a phylogenetic Maximum Likelihood (ML) tree.", "\033[0m")))
    
    
    #Combining OLD database if provided---------------------------------------------------------
    
    if(is.null(old_lib_csv)==FALSE){
      old_lib_file <- read.csv(old_lib_csv)#, row.names=1)
      message(cat(paste0("\n", "\033[0;", 32, "m","Checking old/new library CSV files have same dimensions", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m",paste("No. of columns in old library: ", ncol(old_lib_file), sep=""), "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m",paste("No. of columns in new library: ", ncol(merged.drep1), sep=""), "\033[0m")))
      if(!ncol(old_lib_file)==ncol(merged.drep1)) stop('Library file dimensions are not the same and require manual inspection.', call.=FALSE)
      message(cat(paste0("\033[0;", 32, "m","Same dimensions detected, proceeding...", "\033[0m", "\n")))
      
      #Combine old and new isoLIB files-----------------------
      combined_lib_file <- rbind(old_lib_file, merged.drep1 %>% mutate(representative="")) #old_lib_file must be at top so older seqs take priority
      
      #Check for duplicate names
      combined_lib_file <- combined_lib_file %>% 
        mutate(undup = stringr::str_split_fixed(filename, "_dup", 2)[,1]) %>% 
        mutate(duplicated = duplicated(undup)) %>%
        mutate(duplicated2 = paste(undup, duplicated, sep="___")) %>%
        group_by(duplicated2) %>%
        mutate(filename = ifelse(any(duplicated), paste(undup, "_dup", rep(1:n()), sep=""), filename)) %>%
        ungroup()
      
      if(length((combined_lib_file %>% filter(duplicated==TRUE))$filename) > 0){
        message(cat(paste0("\033[0;", 31, "m", "Duplicated filenames detected: ", paste(unique((combined_lib_file %>% filter(duplicated==TRUE))$undup), collapse="|"), "\033[0m")))
        message(cat(paste0("\033[0;", 30, "m", "Renaming duplicates with suffix '_dup1','_dup2', '_dup3', etc.", "\033[0m")))
        Sys.sleep(2)
        message(cat(paste0("\033[0;", 30, "m", paste((combined_lib_file %>% filter(duplicated==TRUE))$filename, collapse="\n"), "\033[0m")))
      }
      
      combined_lib_file <- combined_lib_file %>% select(-undup, -duplicated, -duplicated2)
      
      #Make unique strain names incase any duplicated-----------------------
      combined_lib_file$unique_filename <- paste(stringr::str_pad(1:nrow(combined_lib_file), 4, pad="0"), combined_lib_file$filename, sep="_")
      
      #New + old sequences
      combined.seqs <- setNames(combined_lib_file$seqs_trim, combined_lib_file$unique_filename)
      #Align sequences
      combined.cluster.aln <- DECIPHER::AlignSeqs(DNAStringSet(combined.seqs))
      #Calculate distance matrix
      combined.cluster.dist <- DECIPHER::DistanceMatrix(combined.cluster.aln, type="matrix")
      #Cluster sequences
      combined.cluster.list <- DECIPHER::TreeLine(myXStringSet = combined.cluster.aln,
                                                  myDistMatrix = combined.cluster.dist,
                                                  cutoff = 1-group_cutoff,
                                                  method = "ML", #UPGMA
                                                  model="JC69",
                                                  maxTime=0.02,
                                                  processors=NULL,
                                                  type="both")
      
      #:::::::::::::::::::::::::::::::::::::::::::::::::::
      #Group sequences based on phylogenetic clustering
      #:::::::::::::::::::::::::::::::::::::::::::::::::::
      if(keep_old_reps==TRUE){
        merged.drep2 <- combined_lib_file %>% 
          mutate(sequence_group.x = combined.cluster.list[[1]]$cluster) %>%
          mutate(row_no = row_number()) %>%
          arrange(Ns_trim) %>%
          arrange(desc(phred_trim)) %>%
          arrange(desc(length_trim)) %>%
          arrange(desc(representative)) %>%
          group_by(sequence_group.x) %>%
          mutate(sequence_group = first(unique_filename)) %>%
          mutate(sequence_group = ifelse(first(representative)=="no", (unique_filename[!grepl("no|yes", representative)])[1], sequence_group)) %>%
          ungroup() %>%
          mutate(sequence_group = sapply(1:length(.$sequence_group),function(x) stringr::str_sub(.$sequence_group[x], 6, nchar(.$sequence_group[x])))) %>%
          mutate(representative = ifelse(sequence_group == filename, "yes", "no")) %>%
          arrange(row_no) %>%
          dplyr::select(sequence_group, date, filename, seqs_trim, phred_trim, Ns_trim, length_trim, closest_match, NCBI_acc, ID,
                        rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, representative)
        
        old.match.comparisons1 <- combined_lib_file %>% filter(representative!="") %>% select(filename, sequence_group) %>% distinct(filename, .keep_all=TRUE)
        old.match.comparisons2 <- merged.drep2 %>% filter(filename %in% tmp$filename) %>% rename("sequence_group2" = "sequence_group") %>% distinct(filename, .keep_all=TRUE) %>% select(filename, sequence_group2)
        
        old.match.comparisons <- cbind(combined_lib_file %>% filter(representative!="") %>% select(filename, sequence_group) %>% distinct(filename, .keep_all=TRUE) %>% {tmp <<- .},
                                       merged.drep2 %>% filter(filename %in% tmp$filename) %>% rename("sequence_group2" = "sequence_group") %>% distinct(filename, .keep_all=TRUE) %>% select(sequence_group2)) %>%
          filter(sequence_group != sequence_group2) #%>% `colnames<-`(c("filename", "Old", "New"))
        
        if(length(old.match.comparisons$filename) > 0){
          message(cat(paste0("\033[0;", 32, "m","Keeping original sequence representatives in place from old library. The following ", length(old.match.comparisons$filename), 
                             " non-representative sequences were assigned to a new group after merging old + new libraries: \n", "\033[0m")))
          message(cat(paste0("\033[0;", 30, "m", "\n",(old.match.comparisons$filename), "\033[0m")))
        } 
      } else {
        merged.drep2 <- combined_lib_file %>% 
          mutate(sequence_group.x = combined.cluster.list[[1]]$cluster) %>%
          mutate(row_no = row_number()) %>%
          arrange(Ns_trim) %>%
          arrange(desc(phred_trim)) %>%
          arrange(desc(length_trim)) %>%
          group_by(sequence_group.x) %>%
          mutate(sequence_group = first(unique_filename)) %>%
          ungroup() %>%
          mutate(sequence_group = sapply(1:length(.$sequence_group),function(x) stringr::str_sub(.$sequence_group[x], 6, nchar(.$sequence_group[x])))) %>%
          mutate(representative = ifelse(sequence_group == filename, "yes", "no")) %>%
          arrange(row_no) %>%
          dplyr::select(sequence_group, date, filename, representative, seqs_trim, phred_trim, Ns_trim, length_trim, closest_match, NCBI_acc, ID,
                        rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species)
        message(cat(paste0("\033[0;", 32, "m","Re-assigning sequence group representatives for old + new libraries.", "\033[0m")))
        
      }
      
      #---Sanity check
      message(cat(paste0("\033[0;", 32, "m","**********Sanity check**********", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m","After merging, the library has a total of ", length(unique(merged.drep2$sequence_group)), 
                         " unique sequence groups detected based on natural clustering on a phylogenetic Maximum Likelihood (ML) tree.", "\033[0m")))
      
    } #-----------end of old library add-in
  } #-----------end of method="maximum_likelihood"
  
  
  ##########################################################################
  ##########################################################################
  ###############                                            ###############
  ###############  START /// method = "closest_species" ///  ###############
  ###############                                            ###############
  ##########################################################################
  ##########################################################################
  
  if(method=="closest_species"){
    #Group sequences based on species-level classification
    merged.drep1 <- new_lib_file %>% 
      mutate(row_no = row_number()) %>% 
      arrange(Ns_trim) %>%
      arrange(desc(phred_trim)) %>%
      arrange(desc(length_trim)) %>%
      arrange(desc(ID)) %>%
      group_by(rank_species) %>%
      mutate(sequence_group = first(unique_filename)) %>%
      ungroup() %>%
      mutate(sequence_group = sapply(1:length(.$sequence_group),function(x) stringr::str_sub(.$sequence_group[x], 6, nchar(.$sequence_group[x])))) %>%
      mutate(representative = ifelse(sequence_group == filename, "yes", "no")) %>%
      arrange(row_no) %>%
      dplyr::select(sequence_group, date, filename, representative, seqs_trim, phred_trim, Ns_trim, length_trim, closest_match, NCBI_acc, ID,
                    rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species)
    
    #---Sanity check
    message(cat(paste0("\033[0;", 32, "m","A total of ", length(unique(merged.drep1$sequence_group)),
                       " unique sequence groups detected based on species-level classification.", "\033[0m")))
    
    
    #Combining OLD database if provided---------------------------------------------------------
    
    if(is.null(old_lib_csv)==FALSE){
      old_lib_file <- read.csv(old_lib_csv)#, row.names=1)
      message(cat(paste0("\n", "\033[0;", 32, "m","Checking old/new library CSV files have same dimensions", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m",paste("No. of columns in old library: ", ncol(old_lib_file), sep=""), "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m",paste("No. of columns in new library: ", ncol(merged.drep1), sep=""), "\033[0m")))
      if(!ncol(old_lib_file)==ncol(merged.drep1)) stop('Library file dimensions are not the same and require manual inspection.', call.=FALSE)
      message(cat(paste0("\033[0;", 32, "m","Same dimensions detected, proceeding...", "\033[0m", "\n")))
      
      #Combine old and new isoLIB files-----------------------
      combined_lib_file <- rbind(old_lib_file, merged.drep1 %>% mutate(representative="")) #old_lib_file must be at top so older seqs take priority
      
      #Check for duplicate names
      combined_lib_file <- combined_lib_file %>% 
        mutate(undup = stringr::str_split_fixed(filename, "_dup", 2)[,1]) %>% 
        mutate(duplicated = duplicated(undup)) %>%
        mutate(duplicated2 = paste(undup, duplicated, sep="___")) %>%
        group_by(duplicated2) %>%
        mutate(filename = ifelse(any(duplicated), paste(undup, "_dup", rep(1:n()), sep=""), filename)) %>%
        ungroup()
      
      if(length((combined_lib_file %>% filter(duplicated==TRUE))$filename) > 0){
        message(cat(paste0("\033[0;", 31, "m", "Duplicated filenames detected: ", paste(unique((combined_lib_file %>% filter(duplicated==TRUE))$undup), collapse="|"), "\033[0m")))
        message(cat(paste0("\033[0;", 30, "m", "Renaming duplicates with suffix '_dup1','_dup2', '_dup3', etc.", "\033[0m")))
        Sys.sleep(2)
        message(cat(paste0("\033[0;", 30, "m", paste((combined_lib_file %>% filter(duplicated==TRUE))$filename, collapse="\n"), "\033[0m")))
      }
      
      combined_lib_file <- combined_lib_file %>% select(-undup, -duplicated, -duplicated2)
      
      #Make unique strain names incase any duplicated-----------------------
      combined_lib_file$unique_filename <- paste(stringr::str_pad(1:nrow(combined_lib_file), 4, pad="0"), combined_lib_file$filename, sep="_")
      
      if(keep_old_reps==TRUE){
        #Group sequences based on species-level classification
        merged.drep2 <- combined_lib_file %>% 
          mutate(row_no = row_number()) %>% 
          arrange(Ns_trim) %>%
          arrange(desc(phred_trim)) %>%
          arrange(desc(length_trim)) %>%
          arrange(desc(ID)) %>%
          arrange(desc(representative)) %>%
          group_by(rank_species) %>%
          mutate(sequence_group = first(unique_filename)) %>%
          ungroup() %>%
          mutate(sequence_group = sapply(1:length(.$sequence_group),function(x) stringr::str_sub(.$sequence_group[x], 6, nchar(.$sequence_group[x])))) %>%
          mutate(representative = ifelse(sequence_group == filename, "yes", "no")) %>%
          arrange(row_no) %>%
          dplyr::select(sequence_group, date, filename, seqs_trim, phred_trim, Ns_trim, length_trim, closest_match, NCBI_acc, ID,
                        rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, representative)
        
        
        old.match.comparisons <- cbind(combined_lib_file %>% filter(representative!="") %>% select(filename, sequence_group) %>% distinct(filename, .keep_all=TRUE) %>% {tmp <<- .},
                                       merged.drep2 %>% filter(filename %in% tmp$filename) %>% rename("sequence_group2" = "sequence_group") %>% distinct(filename, .keep_all=TRUE) %>% select(sequence_group2)) %>% filter(sequence_group != sequence_group2)
        
        if(length(old.match.comparisons$filename) > 0){
          message(cat(paste0("\033[0;", 32, "m","Keeping original sequence representatives in place from old library. The following ", length(old.match.comparisons$filename), 
                             " non-representative sequences were assigned to a new group after merging old + new libraries: \n", "\033[0m")))
          message(cat(paste0("\033[0;", 30, "m", "\n",(old.match.comparisons$filename), "\033[0m")))
        } 
      }else{
        merged.drep2 <- combined_lib_file %>% 
          mutate(row_no = row_number()) %>% 
          arrange(Ns_trim) %>%
          arrange(desc(phred_trim)) %>%
          arrange(desc(length_trim)) %>%
          arrange(desc(ID)) %>%
          group_by(rank_species) %>%
          mutate(sequence_group = first(unique_filename)) %>%
          ungroup() %>%
          mutate(sequence_group = sapply(1:length(.$sequence_group),function(x) stringr::str_sub(.$sequence_group[x], 6, nchar(.$sequence_group[x])))) %>%
          mutate(representative = ifelse(sequence_group == filename, "yes", "no")) %>%
          arrange(row_no) %>%
          dplyr::select(sequence_group, date, filename, seqs_trim, phred_trim, Ns_trim, length_trim, closest_match, NCBI_acc, ID,
                        rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, representative)
        message(cat(paste0("\033[0;", 32, "m","Re-assigning sequence group representatives for old + new libraries.", "\033[0m")))
      }
      
      #---Sanity check
      
      
      message(cat(paste0("\033[0;", 32, "m","**********Sanity check**********", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m","After merging, the library has a total of ", length(unique(merged.drep2$sequence_group)), 
                         " unique sequence groups detected based on species-level classification.", "\033[0m")))
      
    } #-----------end of old library add-in
  } #-----------end of method="closest_species"
  
  
  ##########################################################################
  ##########################################################################
  ###############                                            ###############
  ###############     START /// method = "dark_mode" ///     ###############
  ###############                                            ###############
  ##########################################################################
  ##########################################################################
  if(method=="dark_mode"){
    if(include_warnings == TRUE){
      write.csv(new_lib_file, file=paste(temp_csv, sep=""))
      new_lib_file <- read.csv(paste(temp_csv, sep=""), row.names = 1)
      new_lib_file_path <- paste(temp_csv, sep="")
    } else {
      invisible(capture.output(capture.output(new_lib_file <- tryCatch(
        {
          new_lib_file %>% filter(!grepl("Poor", warning)) %>% filter(!grepl("Fail", decision))
        },
        error = function(e) {
          new_lib_file %>% filter(!grepl("Poor", warning))
        }), type="message"), type="output"))
      if(nrow(new_lib_file)==0)stop("There are no sequences after removing warning sequences. Set 'include_warnings=TRUE' to proceed", call.=FALSE)
      write.csv(new_lib_file, file=paste(temp_csv, sep=""))
      new_lib_file <- read.csv(paste(temp_csv, sep=""), row.names = 1)
      new_lib_file_path <- paste(temp_csv, sep="")
    }
    
    make_fasta((file.path(getwd(), new_lib_file_path)), col_names="unique_filename", col_seqs="seqs_trim", output=paste(temp_fasta, sep="")) # exports as "output.fasta"
    
    #------------------------------------------------- 
    system2(vsearch.path, paste(" --usearch_global ", temp_fasta, " --db ", temp_fasta, " --uc_allhits --uc ", temp_drep, " --id ", group_cutoff, " --query_cov 0.95 --maxaccepts 0 --maxrejects 0 --strand plus ", sep=""), stdout="", stderr="")
    #-------------------------------------------------
    
    #::::::::::::::::::::::::::::
    #Organize search results
    #::::::::::::::::::::::::::::
    drep.results <- read.csv(temp_drep, sep="\t", header = FALSE) %>%
      group_by(V10) %>% mutate(counts = n()) %>% ungroup %>% # Calculate number of matches for each sequence
      .[((sort(.$V8, index.return=TRUE))$ix),] %>% #arrange(desc(as.numeric(V3))) %>%
      arrange(desc(as.numeric(V3))) %>% #V3 is sequence length, sorting from longest to shortest
      arrange(desc(as.numeric(V4))) %>% #V4 is % similarity between query and match, sorting from highest to lowest
      arrange(desc(as.numeric(counts))) %>% #Sorting from sequence with most matches to least matches
      mutate(row_no = row_number()) %>% #distinct(V9, .keep_all=TRUE)
      dplyr::rename("filename" = "V9")
    
    #::::::::::::::::::::::::::::
    # Grouping - Iteration 1
    #::::::::::::::::::::::::::::
    #--------------------------------------------------------------
    unique.groups <- unique(drep.results$V10) #Subtracting list
    match.index <- c()                        #Adding list
    #--------------------------------------------------------------
    while (length(unique.groups) > 0) {
      unique.groups.x <- drep.results %>% filter(V10==unique.groups[1]) %>% filter(!(filename %in% names(match.index)))
      match.index.add <- setNames(unique.groups.x$V10, unique.groups.x$filename) 
      match.index <- c(match.index, match.index.add)
      # check for success
      unique.groups <- unique.groups[!(unique.groups %in% unique.groups[1])]
      unique.groups <- unique.groups[!(unique.groups %in% unique.groups.x$filename)]
    }
    
    #::::::::::::::::::::::::::::
    # Grouping - Iteration 2
    #::::::::::::::::::::::::::::
    #-----------------------------------------------------------------------------
    drep.results.x <- drep.results %>% 
      filter(V10 %in% (unique(paste(match.index)))) %>% #Filtering so that matches in V10 represent only the ref strains  
      #arrange(V10) #%>% #V10 is match column, sorting from oldest to newest added sequences
      arrange(desc(as.numeric(V4))) %>% #V4 is % similarity between query and match, sorting from highest to lowest
      distinct(filename, .keep_all=TRUE)
    
    unlink(file.path(temp_fasta)) #Remove temp FASTA
    unlink(file.path(temp_drep))  #Remove temp DREP
    unlink(file.path(temp_csv))   #Remove temp CSV
    
    #:::::::::::::::::::::::::::::::::::::::::::::
    #Merge results
    #:::::::::::::::::::::::::::::::::::::::::::::
    match.index2 <- setNames(drep.results.x$V10, drep.results.x$filename)
    
    merged.drep <- drep.results %>% 
      mutate(grouping = dplyr::recode(!!!match.index2, filename, .default="")) %>% #Add strain grouping from index
      filter(V10 == grouping) %>%
      merge(.,  new_lib_file, by.x="filename", by.y="unique_filename", all=TRUE) %>% 
      mutate(sequence_group = grouping) %>%
      mutate(representative = ifelse(sequence_group == filename, "yes", "no")) %>%
      mutate(filename = sapply(1:length(.$filename),function(x) stringr::str_sub(.$filename[x], 6, nchar(.$filename[x])))) %>%
      mutate(sequence_group = sapply(1:length(.$sequence_group),function(x) stringr::str_sub(.$sequence_group[x], 6, nchar(.$sequence_group[x]))))
    
    #Subset only columns of interest------------------------------------------------------------------
    
    merged.drep1 <- merged.drep %>% select(sequence_group, date, filename, representative, seqs_trim, phred_trim, Ns_trim, length_trim,
                                           closest_match, NCBI_acc, ID,
                                           rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species)
    
    #---Sanity check
    message(cat(paste0("\033[0;", 32, "m","A total of ", length(unique(merged.drep$grouping)), 
                       " unique sequence groups detected (group_cutoff=", group_cutoff, ")", "\033[0m")))
    
    
    #Combining OLD database if provided---------------------------------------------------------
    if(is.null(old_lib_csv)==FALSE){
      old_lib_file <- read.csv(old_lib_csv)#, row.names=1)
      message(cat(paste0("\n", "\033[0;", 32, "m","Checking old/new library CSV files have same dimensions", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m",paste("No. of columns in old library: ", ncol(old_lib_file), sep=""), "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m",paste("No. of columns in new library: ", ncol(merged.drep1), sep=""), "\033[0m")))
      if(!ncol(old_lib_file)==ncol(merged.drep1)) stop('Library file dimensions are not the same and require manual inspection.', call.=FALSE)
      message(cat(paste0("\033[0;", 32, "m","Same dimensions detected, proceeding...", "\033[0m", "\n")))
      #
      
      #Combine old and new isoLIB files-----------------------
      combined_lib_file <- rbind(old_lib_file, merged.drep1 %>% mutate(representative="")) #old_lib_file must be at top so older seqs take priority
      
      #Check for duplicate names
      combined_lib_file <- combined_lib_file %>% 
        mutate(undup = stringr::str_split_fixed(filename, "_dup", 2)[,1]) %>% 
        mutate(duplicated = duplicated(undup)) %>%
        mutate(duplicated2 = paste(undup, duplicated, sep="___")) %>%
        group_by(duplicated2) %>%
        mutate(filename = ifelse(any(duplicated), paste(undup, "_dup", rep(1:n()), sep=""), filename)) %>%
        ungroup()
      
      if(length((combined_lib_file %>% filter(duplicated==TRUE))$filename) > 0){
        message(cat(paste0("\033[0;", 31, "m", "Duplicated filenames detected: ", paste(unique((combined_lib_file %>% filter(duplicated==TRUE))$undup), collapse="|"), "\033[0m")))
        message(cat(paste0("\033[0;", 30, "m", "Renaming duplicates with suffix '_dup1','_dup2', '_dup3', etc.", "\033[0m")))
        message(cat(paste0("\033[0;", 30, "m", paste((combined_lib_file %>% filter(duplicated==TRUE))$filename, collapse="\n"), "\033[0m")))
      }
      
      combined_lib_file <- combined_lib_file %>% select(-undup, -duplicated, -duplicated2)
      
      #Make unique strain names incase any duplicated-----------------------
      combined_lib_file$unique_filename <- paste(stringr::str_pad(1:nrow(combined_lib_file), 4, pad="0"), combined_lib_file$filename, sep="_")
      
      #Write temp files
      #-------------------------------------------------------
      temp_combined_fasta <- "03_isoLIB_combined_temp.fasta"
      temp_combined_drep <- paste0("03_isoLIB_combined_temp.drep", sep="")
      temp_combined_csv <- "03_isoLIB_combined_temp.csv"
      #-------------------------------------------------------
      write.csv(combined_lib_file, file= temp_combined_csv)
      make_fasta((file.path(getwd(), temp_combined_csv)), col_names="unique_filename", col_seqs="seqs_trim", output=paste(temp_combined_fasta, sep="")) # exports as "output.fasta"
      
      #VSEARCH
      #-------------------------------------------------------
      system2(vsearch.path, paste(" --usearch_global ", temp_combined_fasta, " --db ", temp_combined_fasta, " --uc ", temp_combined_drep, " --id ", group_cutoff, " --query_cov 0.95 --maxaccepts 0 --maxrejects 0 --uc_allhits --strand plus ", sep=""), stdout="", stderr="")
      #-------------------------------------------------------
      
      #::::::::::::::::::::::::::::
      #Organize search results
      #::::::::::::::::::::::::::::
      ref.index <- setNames((combined_lib_file %>% filter(representative!="") %>% filter(filename %in% old_lib_file$filename))$representative,
                            (combined_lib_file %>% filter(representative!="") %>% filter(filename %in% old_lib_file$filename))$unique_filename)
      
      if(keep_old_reps==TRUE){
        drep.results2 <- read.csv(temp_combined_drep, sep="\t", header = FALSE) %>% 
          group_by(V10) %>% mutate(counts = n()) %>% ungroup %>% # Calculate number of matches for each sequence
          .[((sort(.$V8, index.return=TRUE))$ix),] %>% #arrange(desc(as.numeric(V3))) %>%
          arrange(desc(as.numeric(V3))) %>% #V3 is sequence length, sorting from longest to shortest
          arrange(desc(as.numeric(V4))) %>% #V4 is % similarity between query and match, sorting from highest to lowest
          arrange(desc(as.numeric(counts))) %>% #Sorting from sequence with most matches to least matches
          mutate(row_no = row_number()) %>% #distinct(V9, .keep_all=TRUE)
          dplyr::rename("filename" = "V9") %>%
          {drep.results3 <<- .} %>%
          filter(!(V10 %in% names(ref.index)[ref.index=="no"])) #Filtering so that matches in V10 do not contain non-ref strains from 'old_lib_csv' 
      } else {
        drep.results2 <- read.csv(temp_combined_drep, sep="\t", header = FALSE) %>% 
          group_by(V10) %>% mutate(counts = n()) %>% ungroup %>% # Calculate number of matches for each sequence
          .[((sort(.$V8, index.return=TRUE))$ix),] %>% #arrange(desc(as.numeric(V3))) %>%
          arrange(desc(as.numeric(V3))) %>% #V3 is sequence length, sorting from longest to shortest
          arrange(desc(as.numeric(V4))) %>% #V4 is % similarity between query and match, sorting from highest to lowest
          arrange(desc(as.numeric(counts))) %>% #Sorting from sequence with most matches to least matches
          mutate(row_no = row_number()) %>% #distinct(V9, .keep_all=TRUE)
          dplyr::rename("filename" = "V9") %>%
          {drep.results3 <<- .}
      }
      
      #::::::::::::::::::::::::::::
      # Regrouping - Iteration 1
      #::::::::::::::::::::::::::::
      #-----------------------------------------------------------------------------
      unique.groups2 <- sort(unique(drep.results2$V10))   #Subtracting list - #Sorting so that sequences from 'old_lib_csv' are at top and have priority
      match.index2 <- c()                                 #Adding list
      #-----------------------------------------------------------------------------
      while (length(unique.groups2) > 0) {
        match.index2.x <- drep.results2 %>% filter(V10==unique.groups2[1]) %>% filter(!(filename %in% names(match.index2)))
        match.index2.add <- setNames(match.index2.x$V10, match.index2.x$filename) 
        match.index2 <- c(match.index2, match.index2.add)
        # check for success
        unique.groups2 <- unique.groups2[!(unique.groups2 %in% unique.groups2[1])]
        unique.groups2 <- unique.groups2[!(unique.groups2 %in% match.index2.x$filename)]
      }
      #-----------------------------------------------------------------------------
      
      #::::::::::::::::::::::::::::
      # Regrouping - Iteration 2
      #::::::::::::::::::::::::::::
      #-----------------------------------------------------------------------------
      if(keep_old_reps==TRUE){
      drep.results3.x <- drep.results3 %>% 
        filter(V10 %in% unique(paste(match.index2))) %>% #Filtering so that matches in V10 represent only the ref strains  
        arrange(desc(as.numeric(V4))) %>% #V4 is % similarity between query and match, sorting from highest to lowest
        distinct(filename, .keep_all=TRUE) %>% 
        mutate(ref_check1=dplyr::recode(!!!ref.index, filename)) %>%
        mutate(ref_check2=dplyr::recode(!!!ref.index, V10)) %>% 
        mutate(ref_check3= ifelse(ref_check1=="yes" & ref_check2!="yes", filename, "")) %>%
        group_by(V10) %>% 
        mutate(ref_check4 = ifelse(any(ref_check3!=""), (ref_check3[ref_check3!=""])[1], V10)) %>%
        ungroup() %>%
        mutate(V10=ref_check4) %>% dplyr::select(!matches("ref_check"))
      } else {
      drep.results3.x <- drep.results3 %>% 
        filter(V10 %in% unique(paste(match.index2))) %>% #Filtering so that matches in V10 represent only the ref strains  
        arrange(desc(as.numeric(V4))) %>% #V4 is % similarity between query and match, sorting from highest to lowest
        distinct(filename, .keep_all=TRUE)
      }
      
      unlink(file.path(temp_combined_fasta)) #Remove temp FASTA
      unlink(file.path(temp_combined_drep))  #Remove temp DREP
      unlink(file.path(temp_combined_csv))   #Remove temp CSV
      
      #:::::::::::::::::::::::::::::::::::::::::::::
      #Merge results
      #:::::::::::::::::::::::::::::::::::::::::::::
      match.index3 <- setNames(drep.results3.x$V10, drep.results3.x$filename)
      
      combined_lib_file_regroup <- combined_lib_file %>%
        mutate(grouping = dplyr::recode(!!!match.index3, unique_filename, .default="")) %>% #%>% select(-warning) %>% #select(-phylum,-class,-order,-family,-class,-genus, -species_col, -genus_col, -family_col) %>%
        arrange(unique_filename) %>%
        mutate(sequence_group = grouping) %>%
        mutate(representative = ifelse(sequence_group == unique_filename, "yes", "no")) %>%
        #Remove unique identifier characters
        mutate(sequence_group = sapply(1:length(.$sequence_group),function(x) stringr::str_sub(.$sequence_group[x], 6, nchar(.$sequence_group[x]))))
      
      merged.drep2 <- combined_lib_file_regroup %>% select(sequence_group, date, filename, representative, seqs_trim, phred_trim, Ns_trim, length_trim,
                                                           closest_match, NCBI_acc, ID,
                                                           rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species)
      
      #---Sanity check
      old.match.comparisons <- cbind(combined_lib_file %>% filter(representative!="") %>% select(filename, sequence_group) %>% distinct(filename, .keep_all=TRUE) %>% {tmp <<- .},
                                     merged.drep2 %>% filter(filename %in% tmp$filename) %>% rename("sequence_group2" = "sequence_group") %>% distinct(filename, .keep_all=TRUE) %>% select(sequence_group2)) %>%
        filter(sequence_group != sequence_group2) #%>% `colnames<-`(c("filename", "Old", "New"))
      
      message(cat(paste0("\033[0;", 32, "m","**********Sanity check**********", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m","After merging, the library has a total of ", length(unique(merged.drep2$sequence_group)), 
                         " unique sequence groups (group_cutoff=", group_cutoff, ")", "\033[0m")))
      
      if(length(old.match.comparisons$filename) > 0){
        message(cat(paste0("\033[0;", 32, "m","The following ", length(old.match.comparisons$filename), 
                           " strains were assigned to new reference sequence groups after merging old + new libraries: \n", "\033[0m")))
        message(cat(paste0("\033[0;", 30, "m", "\n",(old.match.comparisons$filename), "\033[0m")))
      }
      
    } #-----------end of old library add-in
  } #-----------end of method="dark_mode"
  
  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################
  
  #::::::::::::::::::::::::::::::::::::::::::::::
  # Set cutoffs for taxonomic rank demarcation
  #::::::::::::::::::::::::::::::::::::::::::::::
  
  
  merged.drep1 <- merged.drep1 %>%
    mutate(phylum_threshold= phylum_threshold) %>%
    mutate(class_threshold= class_threshold) %>%
    mutate(order_threshold= order_threshold) %>%
    mutate(family_threshold= family_threshold) %>%
    mutate(genus_threshold= genus_threshold) %>%
    mutate(species_threshold= species_threshold)
  
  if(is.null(old_lib_csv)==FALSE){
    merged.drep2 <- merged.drep2 %>%
      mutate(phylum_threshold= phylum_threshold) %>%
      mutate(class_threshold= class_threshold) %>%
      mutate(order_threshold= order_threshold) %>%
      mutate(family_threshold= family_threshold) %>%
      mutate(genus_threshold= genus_threshold) %>%
      mutate(species_threshold= species_threshold)
  }
  
  ###############
  # Output
  ###############
  
  #Set S4 class for outputs
  isolib.S4.1 <- df_to_isoLIB(merged.drep1)
  if(is.null(old_lib_csv)==FALSE){isolib.S4.2 <- df_to_isoLIB(merged.drep2)}
  
  #--------------------------------------------------
  #Export messages
  message(cat(paste0("\033[97;", 40, "m","'isoLIB' steps completed. Exporting files...", "\033[0m")))
  
  #Export HTML file----------------------------------------------------------------------
  if(export_html==TRUE){
    if(is.null(old_lib_csv)==TRUE){
      output_name <- "03_isoLIB_results.html"
      export_html(isolib.S4.1, method = method, group_cutoff=group_cutoff)
      #merged.df3_react <- export_html(out_df, path, output=output_name)
      message(cat(paste0("\033[97;", 40, "m","HTML results exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m", "\n")))
    }
    if(is.null(old_lib_csv)==FALSE){
      output_name <- "03_isoLIB_results_merged.html"
      export_html(isolib.S4.2, method = method, group_cutoff=group_cutoff)
      #merged.df3_react <- export_html(out_df, path, output=output_name)
      message(cat(paste0("\033[97;", 40, "m","HTML results exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m", "\n")))
    }
  } else if(export_html==FALSE){}
  
  #export CSV files----------------------------------------------------------------------
  if(export_csv==TRUE){
    if(is.null(old_lib_csv)==TRUE){
      output_name <- "03_isoLIB_results.csv"
      csv_output <- S4_to_dataframe(isolib.S4.1) %>% select(sequence_group, date, filename, representative, seqs_trim, phred_trim, Ns_trim, length_trim,
                                                            closest_match, NCBI_acc, ID,
                                                            rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species)
      write.csv(csv_output, file=paste(output_name, sep=""), row.names = FALSE)
      message(cat(paste0("\033[97;", 40, "m","CSV file exported:", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m",
                         "\033[0;", 31, "m", "  <--- Final strain library output","\033[0m")))
    }
    if(is.null(old_lib_csv)==FALSE){
      output_name <- "03_isoLIB_results_merged.csv"
      csv_output <- S4_to_dataframe(isolib.S4.2) %>% select(sequence_group, date, filename, representative, seqs_trim, phred_trim, Ns_trim, length_trim,
                                                            closest_match, NCBI_acc, ID,
                                                            rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species)
      write.csv(csv_output, file=paste(output_name, sep=""), row.names = FALSE)
      message(cat(paste0("\033[97;", 40, "m","CSV file exported:", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m",
                         "\033[0;", 31, "m", "  <--- Final MERGED strain library output","\033[0m")))
    }
  } else if(export_csv==FALSE){}
  
  if(is.null(old_lib_csv)==TRUE){return(isolib.S4.1)}
  if(is.null(old_lib_csv)==FALSE){return(isolib.S4.2)}
  
} #end of function

###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################
###########################################################################################################





