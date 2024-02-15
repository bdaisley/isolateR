#' @title Generate new strain library or add to existing one.
#' @name isoLIB
#' @rdname isoLIB
#' @description This function creates a strain library by grouping closely related strains of interest based on sequence similarity.
#' For adding new sequences to an already-established strain library, specify the .CSV file path of the older strain library using the 'old_lib_csv" parameter.
#' @export
#' @param input Path of CSV output file from isoTAX step.
#' @param old_lib_csv Optional: Path of CSV output isoLIB file or combined isoLIB file from previous run(s)
#' @param export_html (Default=TRUE) Output the results as an HTML file
#' @param export_csv (Default=TRUE) Output the results as a CSV file.
#' @param include_warnings (Default=FALSE) Whether or not to keep sequences with poor alignment warnings from Step 2 'isoTAX' function. Set TRUE to keep warning sequences, and FALSE to remove warning sequences.
#' @param strain_group_cutoff (Default=0.995) Similarity cutoff (0-1) for delineating between strain groups. 1 = 100% identical/0.995=0.5% difference/0.95=5.0% difference/etc.
#' @param phylum_cutoff Percent cutoff for phylum rank demarcation
#' @param class_cutoff Percent cutoff for class rank demarcation
#' @param order_cutoff Percent cutoff for order rank demarcation
#' @param family_cutoff Percent cutoff for family rank demarcation
#' @param genus_cutoff Percent cutoff for genus rank demarcation
#' @param species_cutoff Percent cutoff for species rank demarcation
#' @seealso \code{\link{isoTAX}}, \code{\link{isoLIB}}
#' @return Returns an isoLIB class object. Default taxonomic cutoffs for phylum (75.0), class (78.5), order (82.0), family (86.5), genus (96.5), and species (98.7) demarcation are based on Yarza et al. 2014, Nature Reviews Microbiology (DOI:10.1038/nrmicro3330)
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
                   export_html=TRUE,
                   export_csv=TRUE,
                   include_warnings=FALSE,
                   strain_group_cutoff = 0.995,
                   phylum_cutoff=75.0,
                   class_cutoff=78.5,
                   order_cutoff=82.0,
                   family_cutoff=86.5,
                   genus_cutoff=96.5,
                   species_cutoff=98.7){

  #Input checks------------------------------------------------------------

  if(!is.null(strain_group_cutoff)){
    if(!is.numeric(strain_group_cutoff)) stop("'strain_group_cutoff' is not numeric. Set between 0-1 (1=no difference i.e., identical sequences)", call.=FALSE)
    if(strain_group_cutoff <0 | strain_group_cutoff >1) stop("Wrong 'strain_group_cutoff' format. Set between 0-1 (1=no difference i.e., identical sequences)", call.=FALSE)
  }

  #Set paths------------------------------------------------------------

	new_lib <- stringr::str_replace_all(input, '\\\\', '/')
	new_lib_path <- paste(unlist(strsplit(new_lib, '/'))[1:(length(unlist(strsplit(new_lib, '/')))-1)],collapse="/")
	new_lib_file <- read.csv(new_lib) #, row.names = 1) removing row.names setting
	setwd(new_lib_path)
	
  #-------------------------------------------------
  #:::::::::::::::::::::::::::
  #Download VSEARCH software
  #:::::::::::::::::::::::::::

  path = new_lib_path
  vsearch.path.dl <- file.path(system.file("", package="isolateR"), "vsearch")
  suppressWarnings(dir.create(vsearch.path.dl))
  vsearch_files <- stringr::str_subset(dir(vsearch.path.dl, full.names = FALSE), 'vsearch')


  if(paste(isolateR::get_os())=="windows"){
    if(identical(vsearch_files, character(0))){
      message(cat(paste0("\n", "\033[0;", 32, "m","Operating system is ---> Windows-based <---", "\033[0m")))
      download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-win-x86_64.zip", file.path(vsearch.path.dl, 'vsearch-2.23.0-win-x86_64.zip'), mode='wb')
      unzip(file.path(vsearch.path.dl,"vsearch-2.23.0-win-x86_64.zip"),  exdir=file.path(vsearch.path.dl))
      file.copy(file.path(vsearch.path.dl, "vsearch-2.23.0-win-x86_64/bin/vsearch.exe"), file.path(vsearch.path.dl, "vsearch-2.23.0.exe"), overwrite=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-win-x86_64"),recursive=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-win-x86_64.zip"),recursive=TRUE)
      message(cat(paste0("\n", "\033[0;", 32, "m","Download complete. VSEARCH 2.23.0 has been installed.", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0.exe")
    } else {
      message(cat(paste0("\n", "\033[0;", 32, "m","Operating system is ---> Windows-based <---", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m","VSEARCH already downloaded", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0.exe")
    }
  }

  if(paste(isolateR::get_os())=="osx-mac"){
    if(identical(vsearch_files, character(0))){
      message(cat(paste0("\033[0;", 32, "m","Operating system is ---> MacOS-based <---", "\033[0m")))
      download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-macos-universal.tar.gz", file.path(vsearch.path.dl, 'vsearch-2.23.0-macos-universal.tar.gz'), mode='wb')
      untar(file.path(vsearch.path.dl,"vsearch-2.23.0-macos-universal.tar.gz"),  exdir=file.path(vsearch.path.dl))
      file.copy(file.path(vsearch.path.dl, "vsearch-2.23.0-macos-universal/bin/vsearch"), file.path(vsearch.path.dl, "vsearch-2.23.0_macos"), overwrite=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-macos-universal"),recursive=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-macos-universal.tar.gz"),recursive=TRUE)
      message(cat(paste0("\033[0;", 32, "m","Download complete. VSEARCH 2.23.0 has been installed.", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0_macos")
    } else {
      message(cat(paste0("\033[0;", 32, "m","Operating system is ---> MacOS-based <---", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m","VSEARCH already downloaded.", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0_macos")
    }
  }

  if(paste(isolateR::get_os())=="linux"){
    if(identical(vsearch_files, character(0))){
      message(cat(paste0("\033[0;", 32, "m","Operating system is ---> Linux-based <---", "\033[0m")))
      download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-linux-x86_64.tar.gz", file.path(vsearch.path.dl, 'vsearch-2.23.0-linux-x86_64.tar.gz'), mode='wb')
      untar(file.path(vsearch.path.dl, "vsearch-2.23.0-linux-x86_64.tar.gz"),  exdir=file.path(vsearch.path.dl))
      file.copy(file.path(vsearch.path.dl, "vsearch-2.23.0-linux-x86_64/bin/vsearch"), file.path(vsearch.path.dl, "vsearch-2.23.0"), overwrite=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-linux-x86_64"),recursive=TRUE)
      unlink(file.path(vsearch.path.dl,"vsearch-2.23.0-linux-x86_64.tar.gz"),recursive=TRUE)
      message(cat(paste0("\033[0;", 32, "m","Download complete. VSEARCH 2.23.0 has been installed.", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0")
    } else {
      message(cat(paste0("\033[0;", 32, "m","Operating system is ---> Linux-based <---", "\033[0m")))
      message(cat(paste0("\033[0;", 32, "m","VSEARCH already downloaded.", "\033[0m", "\n")))
      vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0")
    }

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
  system2(vsearch.path, paste(" --usearch_global ", temp_fasta, " --db ", temp_fasta, " --uc_allhits --uc ", temp_drep, " --id ", strain_group_cutoff, " --query_cov 0.95 --maxaccepts 0 --maxrejects 0 --strand plus ", sep=""), stdout="", stderr="")
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
    # do something
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
  
  #-----------------------------------------------------------------------------
  unique.groups2 <- sort(unique(drep.results.x$V10))    #Subtracting list
  match.index2 <- c()                                   #Adding list
  #-----------------------------------------------------------------------------
  while (length(unique.groups2) > 0) {
    # do something
    unique.groups2.x <- drep.results.x %>% filter(V10==unique.groups2[1]) %>% filter(!(filename %in% names(match.index2)))
    match.index2.add <- setNames(unique.groups2.x$V10, unique.groups2.x$filename) 
    match.index2 <- c(match.index2, match.index2.add)
    # check for success
    unique.groups2 <- unique.groups2[!(unique.groups2 %in% unique.groups2[1])]
    unique.groups2 <- unique.groups2[!(unique.groups2 %in% unique.groups2.x$filename)]
  }
  #-----------------------------------------------------------------------------

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
    mutate(strain_group = grouping) %>%
    mutate(ref_strain = ifelse(strain_group == filename, "yes", "no")) %>%
    mutate(filename = sapply(1:length(.$filename),function(x) stringr::str_sub(.$filename[x], 6, nchar(.$filename[x])))) %>%
    mutate(strain_group = sapply(1:length(.$strain_group),function(x) stringr::str_sub(.$strain_group[x], 6, nchar(.$strain_group[x]))))
  
  #Subset only columns of interest------------------------------------------------------------------

  merged.drep1 <- merged.drep %>% select(strain_group, date, filename, seqs_trim, phred_trim, Ns_trim, length_trim,
                                              closest_match, NCBI_acc, ID,
                                              rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, ref_strain)

  #:::::::::::::::::::::::::::::::::::::::::::::
  #Sanity check
  #:::::::::::::::::::::::::::::::::::::::::::::
  
  message(cat(paste0("\033[0;", 32, "m","A total of ", length(unique(merged.drep$grouping)), 
                     " unique sequence groups detected (strain_group_cutoff=", strain_group_cutoff, ")", "\033[0m")))
  
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #Combining OLD database if provided---------------------------------------------------------
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

  if(is.null(old_lib_csv)==FALSE){
    old_lib_file <- read.csv(old_lib_csv)#, row.names=1)
    message(cat(paste0("\n", "\033[0;", 32, "m","Checking old/new library CSV files have same dimensions", "\033[0m")))
    message(cat(paste0("\033[0;", 32, "m",paste("No. of columns in old library: ", ncol(old_lib_file), sep=""), "\033[0m")))
    message(cat(paste0("\033[0;", 32, "m",paste("No. of columns in new library: ", ncol(merged.drep1), sep=""), "\033[0m")))
    if(!ncol(old_lib_file)==ncol(merged.drep1)) stop('Library file dimensions are not the same and require manual inspection.', call.=FALSE)
    message(cat(paste0("\033[0;", 32, "m","Same dimensions detected, proceeding...", "\033[0m", "\n")))
    #

    #Combine old and new isoLIB files-----------------------
    combined_lib_file <- rbind(old_lib_file, merged.drep1 %>% mutate(ref_strain="")) #old_lib_file must be at top so older seqs take priority
    
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
    system2(vsearch.path, paste(" --usearch_global ", temp_combined_fasta, " --db ", temp_combined_fasta, " --uc ", temp_combined_drep, " --id ", strain_group_cutoff, " --query_cov 0.95 --maxaccepts 0 --maxrejects 0 --uc_allhits --strand plus ", sep=""), stdout="", stderr="")
    #-------------------------------------------------------
    
    #::::::::::::::::::::::::::::
    #Organize search results
    #::::::::::::::::::::::::::::
    ref.index <- setNames((combined_lib_file %>% filter(ref_strain!="") %>% filter(filename %in% old_lib_file$filename))$ref_strain,
                          (combined_lib_file %>% filter(ref_strain!="") %>% filter(filename %in% old_lib_file$filename))$unique_filename)
    
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

    #::::::::::::::::::::::::::::
    # Regrouping - Iteration 1
    #::::::::::::::::::::::::::::
    #-----------------------------------------------------------------------------
    unique.groups2 <- sort(unique(drep.results2$V10))   #Subtracting list - #Sorting so that sequences from 'old_lib_csv' are at top and have priority
    match.index2 <- c()                                 #Adding list
    #-----------------------------------------------------------------------------
    while (length(unique.groups2) > 0) {
      # do something
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
    drep.results3.x <- drep.results3 %>% 
      #filter(!(filename %in% ref.index.new)) %>%
      filter(V10 %in% unique(paste(match.index2))) %>% #Filtering so that matches in V10 represent only the ref strains  
      #arrange(V10) %>% #V10 is match column, sorting from oldest to newest added sequences
      arrange(desc(as.numeric(V4))) %>% #V4 is % similarity between query and match, sorting from highest to lowest
      distinct(filename, .keep_all=TRUE)

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
      mutate(strain_group = grouping) %>%
      mutate(ref_strain = ifelse(strain_group == unique_filename, "yes", "no")) %>%
      #Remove unique identifier characters
      mutate(strain_group = sapply(1:length(.$strain_group),function(x) stringr::str_sub(.$strain_group[x], 6, nchar(.$strain_group[x]))))
    
    merged.drep2 <- combined_lib_file_regroup %>% select(strain_group, date, filename, seqs_trim, phred_trim, Ns_trim, length_trim,
                                                closest_match, NCBI_acc, ID,
                                                rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, ref_strain)
    
    #:::::::::::::::::::::::::::::::::::::::::::::
    #Sanity check
    #:::::::::::::::::::::::::::::::::::::::::::::
    
    old.match.comparisons <- cbind(combined_lib_file %>% filter(ref_strain!="") %>% select(filename, strain_group) %>% distinct(filename, .keep_all=TRUE) %>% {tmp <<- .},
                                   merged.drep2 %>% filter(filename %in% tmp$filename) %>% rename("strain_group2" = "strain_group") %>% distinct(filename, .keep_all=TRUE) %>% select(strain_group2)) %>%
      filter(strain_group != strain_group2) #%>% `colnames<-`(c("filename", "Old", "New"))
    
    message(cat(paste0("\033[0;", 32, "m","**********Sanity check**********", "\033[0m")))
    message(cat(paste0("\033[0;", 32, "m","After merging, the library has a total of ", length(unique(merged.drep2$strain_group)), 
                       " unique sequence groups (strain_group_cutoff=", strain_group_cutoff, ")", "\033[0m")))
    
    if(length(old.match.comparisons$filename) > 0){
      message(cat(paste0("\033[0;", 32, "m","The following ", length(old.match.comparisons$filename), 
                         " strains were assigned to new reference sequence groups after merging: \n", "\033[0m")))
      message(cat(paste0("\033[0;", 30, "m", "\n",(old.match.comparisons$filename), "\033[0m")))
    }
    
    } #-----------end of old library add-in
    
    
  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################
  ###########################################################################################################

  #::::::::::::::::::::::::::::::::::::::::::::::
  # Set cutoffs for taxonomic rank demarcation
  #::::::::::::::::::::::::::::::::::::::::::::::
  
  
  merged.drep1 <- merged.drep1 %>%
    mutate(phylum_cutoff= phylum_cutoff) %>%
    mutate(class_cutoff= class_cutoff) %>%
    mutate(order_cutoff= order_cutoff) %>%
    mutate(family_cutoff= family_cutoff) %>%
    mutate(genus_cutoff= genus_cutoff) %>%
    mutate(species_cutoff= species_cutoff)

  if(is.null(old_lib_csv)==FALSE){
    merged.drep2 <- merged.drep2 %>%
      mutate(phylum_cutoff= phylum_cutoff) %>%
      mutate(class_cutoff= class_cutoff) %>%
      mutate(order_cutoff= order_cutoff) %>%
      mutate(family_cutoff= family_cutoff) %>%
      mutate(genus_cutoff= genus_cutoff) %>%
      mutate(species_cutoff= species_cutoff)
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
      export_html(isolib.S4.1)
      #merged.df3_react <- export_html(out_df, path, output=output_name)
      message(cat(paste0("\033[97;", 40, "m","HTML results exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m", "\n")))
    }
    if(is.null(old_lib_csv)==FALSE){
      output_name <- "03_isoLIB_results_merged.html"
      export_html(isolib.S4.2)
      #merged.df3_react <- export_html(out_df, path, output=output_name)
      message(cat(paste0("\033[97;", 40, "m","HTML results exported: ", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m", "\n")))
    }
    } else if(export_html==FALSE){}

  #export CSV files----------------------------------------------------------------------
  if(export_csv==TRUE){
    if(is.null(old_lib_csv)==TRUE){
      output_name <- "03_isoLIB_results.csv"
      csv_output <- S4_to_dataframe(isolib.S4.1) %>% select(strain_group, date, filename,  seqs_trim, phred_trim, Ns_trim, length_trim,
                                                            closest_match, NCBI_acc, ID,
                                                            rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, ref_strain)
      write.csv(csv_output, file=paste(output_name, sep=""), row.names = FALSE)
      message(cat(paste0("\033[97;", 40, "m","CSV file exported:", "\033[0m",
                         "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m",
                         "\033[0;", 31, "m", "  <--- Final strain library output","\033[0m")))
    }
    if(is.null(old_lib_csv)==FALSE){
      output_name <- "03_isoLIB_results_merged.csv"
      csv_output <- S4_to_dataframe(isolib.S4.2) %>% select(strain_group, date, filename,  seqs_trim, phred_trim, Ns_trim, length_trim,
                                                          closest_match, NCBI_acc, ID,
                                                          rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, ref_strain)
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





