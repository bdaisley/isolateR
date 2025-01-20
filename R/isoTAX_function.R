#' @name isoTAX
#' @title Classify taxonomy of sequences after quality trimming steps.
#' @rdname isoTAX
#' @description  This function performs taxonomic classification steps by searching query Sanger sequences against specified database of interest. Takes CSV input files, extracts FASTA-formatted query sequences and performs global alignment against specified database of interest via Needleman-Wunsch algorithm by wrapping the --usearch_global command implemented in VSEARCH. Default taxonomic rank cutoffs for 16S rRNA gene sequences are based on Yarza et al. 2014, Nat Rev Microbiol.
#' @export
#' @param input Path of either 1) CSV output file from isoQC step, or 2) a FASTA formatted file.
#' If input is a FASTA file, the sequence(s) will be converted and saved as an isoQC-formatted output file in the current 
#' working directory ("isolateR_output/01_isoQC_mock_table.csv"). Sequence date, name, length, and number of ambiguous bases (Ns) 
#' will be calculated from the input file and used to populate the relevant columns. Phred quality scores (phred_trim) will be 
#' set to the maximum value (60) and the remaining columns will be populated with mock data to allow compatibility with 
#' the isoTAX function. The main purpose of this output file is for flexibility and to allow users to edit/modify the sequence 
#' metadata before continuing with subsequent steps.
#' @param export_html (Default=TRUE) Output the results as an HTML file
#' @param export_csv (Default=TRUE) Output the results as a CSV file.
#' @param export_blast_table (Default=FALSE) Output the results as a tab-separated BLAST-like hits table.
#' @param quick_search (Default=FALSE) Whether or not to perform a comprehensive database search (i.e. optimal global alignment).
#' If TRUE, performs quick search equivalent to setting VSEARCH parameters "--maxaccepts 100 --maxrejects 100".
#' If FALSE, performs comprehensive search equivalent to setting VSEARCH parameters "--maxaccepts 0 --maxrejects 0"
#' @param db (Default="16S_bac") Select database option(s) including "16S" (for searching against the NCBI Refseq targeted loci 16S rRNA database),
#' "ITS" (for searching against the NCBI Refseq targeted loci ITS  database. For combined databases in cases where input sequences are derived from
#' bacteria and fungi, select "16S|ITS". Setting to anything other than db=NULL or db="custom" causes 'db.path' parameter to be ignored.
#' @param db_path Path of FASTA-formatted database sequence file. Ignored if 'db' parameter is set to anything other than NULL or "custom".
#' @param vsearch_path Path of VSEARCH software if manually downloaded in a custom directory. If NULL (Default), will attempt automatic download.
#' @param iddef Set pairwise identity definition as per VSEARCH definitions (Default=2, and is recommended for highest taxonomic accuracy)
#' (0) CD-HIT definition: (matching columns) / (shortest sequence length).
#' (1) Edit distance: (matching columns) / (alignment length).
#' (2) Edit distance excluding terminal gaps (default definition).
#' (3) Marine Biological Lab definition counting each gap opening (internal or terminal) as a single mismatch, whether or not the gap was extended: 1.0- ((mismatches + gap openings)/(longest sequence length)).
#' (4) BLAST definition, equivalent to --iddef 1 for global pairwise alignments.
#' @param phylum_threshold Percent sequence similarity threshold for phylum rank demarcation
#' @param class_threshold Percent sequence similarity threshold for class rank demarcation
#' @param order_threshold Percent sequence similarity threshold for order rank demarcation
#' @param family_threshold Percent sequence similarity threshold for family rank demarcation
#' @param genus_threshold Percent sequence similarity threshold for genus rank demarcation
#' @param species_threshold Percent sequence similarity threshold for species rank demarcation
#' @seealso \code{\link{isoQC}}, \code{\link{isoLIB}}, \code{\link{search_db}}
#' @return Returns taxonomic classification table of class isoTAX. Default taxonomic cutoffs for phylum (75.0), class (78.5), order (82.0), family (86.5), genus (94.5), and species (98.7) demarcation are based on Yarza et al. 2014, Nature Reviews Microbiology (DOI:10.1038/nrmicro3330)
#' @importFrom stringr str_subset
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_split_fixed
#' @importFrom stringr str_remove
#' @importFrom svMisc progress
#' @importFrom xmlconvert xml_to_df
#' @importFrom rentrez entrez_fetch
#' @importFrom utils download.file
#' @importFrom utils untar
#' @importFrom R.utils gunzip
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

#' #Show summary statistics
#' isoTAX.S4

isoTAX <- function(input=NULL,
                   export_html=TRUE,
                   export_csv=TRUE,
                   export_blast_table=FALSE,
                   quick_search=FALSE,
                   db="16S_bac",
                   db_path=NULL,
                   vsearch_path=NULL,
                   iddef=2,
                   phylum_threshold=75.0,
                   class_threshold=78.5,
                   order_threshold=82.0,
                   family_threshold=86.5,
                   genus_threshold=94.5,
                   species_threshold=98.7
                   ){


  #Setting file paths-------------------------------------------------
  input <- stringr::str_replace_all(input, '\\\\', '/')
  path <- paste(unlist(strsplit(input, '/'))[1:(length(unlist(strsplit(input, '/')))-1)], collapse="/")
  setwd(path)
  
  #:::::::::::::::::::::::::::::::::::::::::
  #Evaluate input file type
  #:::::::::::::::::::::::::::::::::::::::::
  
  #Check if input is file (proceed) or directory (stop)
  suppressWarnings(invisible(
    capture.output(
      capture.output(
        file.check <- try(strsplit(readLines(input), "(?=.)", perl=TRUE)[[1]][1]), 
        type="message"), 
      type="output")))
  
  if(is(file.check, 'try-error')) stop('Input requires a file, not a directory. Please specify either an isoQC output file or FASTA file to proceed.', call.=FALSE)

  #If FASTA format
  if(file.check==">"){
    message(cat(paste0("\033[47;", 32, "m","FASTA formatted file detected as input. Converting...", "\033[0m")))
    message(cat(paste0("\033[47;", 32, "m","Note: Phred quality values will be set to the maximum (60) when using FASTA files as input.", "\033[0m")))
    
    fasta.in <- Biostrings::readBStringSet(input)
    
    fasta.date <- stringr::str_split_fixed((file.info(input))$ctime, " ", 2)[,1]
    fasta.date <- paste(stringr::str_pad(stringr::str_split_fixed(fasta.date, "-", 3)[,1], 2, pad = "0"),
                        stringr::str_pad(stringr::str_split_fixed(fasta.date, "-", 3)[,2], 2, pad = "0"),
                        stringr::str_pad(stringr::str_split_fixed(fasta.date, "-", 3)[,3], 2, pad = "0"), sep="_")
    
    entry.fasta <- data.frame(date=fasta.date,
                              filename=gsub(" ", "_",names(fasta.in)),
                              seqs_raw = rep("NNN", length(fasta.in)),
                              phred_raw = rep(0, length(fasta.in)),
                              Ns_raw = rep(0, length(fasta.in)),
                              length_raw = rep(0, length(fasta.in)),
                              phred_spark_raw = paste(rep(0.1, length(fasta.in)), collapse="_"),
                              seqs_trim = paste(fasta.in),
                              phred_trim = rep(60, length(fasta.in)),
                              Ns_trim = stringr::str_count(paste(fasta.in), "[Nn]"),
                              length_trim = nchar(paste(fasta.in)),
                              decision = rep("Pass", length(fasta.in)) )
        
    raw.input <- entry.fasta
    
    #Add isolateR_output directory if it doesn't exist and adjust path accordingly
    if(paste(unlist(strsplit(input, '/'))[(length(unlist(strsplit(input, '/')))-1)], collapse="/") != "isolateR_output"){
      path <- file.path(path, "isolateR_output")
      suppressWarnings(dir.create(path))
      suppressWarnings(dir.create(file.path(path, "lib")))
      setwd(path)
    }
    write.csv(raw.input, file=file.path(path, "01_isoQC_mock_table.csv"), row.names = FALSE)
  }
  
  #Check if isoQC format (proceed) or if file format not recognized (stop)
  if(file.check!=">"){
    entry.csv <- try(read.csv(input))
    if(is(entry.csv, 'try-error')) stop('Input file not recognized. Should be either isoQC output or FASTA format file.', call.=FALSE)
    raw.input <- entry.csv
  }
  
  #:::::::::::::::::::::::::::::::::::::::::
  #Filtering input file, catch duplicates
  #:::::::::::::::::::::::::::::::::::::::::

  #raw.input <- read.csv(input)
    
  #Make unique strain names incase any duplicated-----------------------
    
  raw.input <- raw.input %>% arrange(desc(length_trim))
    
  #Check for duplicate names, and add _dup# for duplicates and remove spaces from filename
  
  raw.input <- raw.input %>%
    mutate(filename = sub(" ", "_", filename)) %>%
    mutate(filename = stringr::str_remove(filename,'_[:alpha:][:digit:][:digit:]\\.ab1$')) %>%
    group_by(filename) %>%
    mutate(dup_count = row_number()) %>%
    mutate(filename = ifelse(dup_count > 1, paste(filename, "_dup", dup_count - 1, sep = ""), filename)) %>%
    mutate(duplicated = dup_count > 1) %>%
    select(-dup_count) %>%
    ungroup()
  
  if(length((raw.input %>% filter(duplicated==TRUE))$filename) > 0){
    message(cat(paste0("\033[0;", 31, "m", "Duplicated filenames detected: ", paste(unique((raw.input %>% filter(duplicated==TRUE))$filename), collapse="|"), "\033[0m")))
    message(cat(paste0("\033[0;", 30, "m", "Renaming duplicates with suffix '_dup1','_dup2', '_dup3', etc.", "\033[0m")))
    message(cat(paste0("\033[0;", 30, "m",  paste((raw.input %>% filter(duplicated==TRUE))$filename, collapse="\n"), "\033[0m")))
  }
  
  isoTAX.input <- raw.input %>% select(-duplicated)
  
  #:::::::::::::::::::::::::::::::::::::::::
  #Stage paths for taxonomic classification
  #:::::::::::::::::::::::::::::::::::::::::
  message(cat(paste0("\033[0;", 32, "m","Staging files", "\033[0m")))
  
  #Stage query FASTA file
  query.fasta.path <- "02_isoTAX_query.fasta"
  df <- as.data.frame(isoTAX.input)
  col_names <- "filename"
  col_seqs <- "seqs_trim"
  output <- query.fasta.path
  
  #removed call to make_fasta, running same script here so updated df from above is used as input
  ss <- Biostrings::DNAStringSet(stats::setNames(df[,paste(col_seqs)], df[,paste(col_names)]))
  ss <- Biostrings::DNAStringSet(c(df[,col_seqs]))
  names(ss) <- paste(df[,col_names])
  fasta.output <- file.path(path, output)
  Biostrings::writeXStringSet(ss, fasta.output)
  
  
  #Output - UC table
  uc.out <- "02_isoTAX_output.uc"
  #Output - BLAST6 table
  b6.out <- "02_isoTAX_output.b6o"
  
  uc.results <- search_db(query.path = query.fasta.path,
                          uc.out = uc.out,
                          b6.out = b6.out,
                          path = path,
                          quick_search = quick_search,
                          db = db,
                          db_path = db_path,
                          vsearch_path= vsearch_path,
                          iddef=iddef,
                          keep_temp_files=export_blast_table)
  
  if(export_blast_table==TRUE){
    blast.table <- read.csv("temp_vsearch/02_isoTAX_output.b6o", sep="\t", header = FALSE)
    unlink("temp_vsearch",recursive=TRUE)
  }
  
  message(cat(paste0("\033[97;", 40, "m","Organizing results", "\033[0m")))
  
  #-------------------------------------------------
  
  query.info <- isoTAX.input
  query.info.list1 <- setNames(query.info$phred_trim, query.info$filename)
  query.info.list2 <- setNames(query.info$date, query.info$filename)
  
  message(cat(paste0("\033[97;", 40, "m","Merging results...", "\033[0m")))
  
  merged.df <- merge(query.info, uc.results, by.x="filename", by.y="V9", all=TRUE) %>%
    mutate(V10 = gsub(" ", ";", V10)) %>%
    mutate(NCBI_acc = stringr::str_split_fixed(.$V10, ";", 8)[,1]) %>%
    dplyr::rename("ID" = "V4", "closest_match" = "V10") %>%
    mutate(NCBI_acc = ifelse(length_trim < 20, "No_match", NCBI_acc)) %>%
    mutate(closest_match = ifelse(length_trim < 20, "No_match", closest_match)) %>%
    mutate(ID = ifelse(length_trim < 20, 51, ID))
  
  message(cat(paste0("\033[0;", 32, "m","Done.", "\033[0m")))
  
  #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  if(db!="custom"){
    #------------------------------------------------- library(rentrez)
    message(cat(paste0("\033[97;", 40, "m","Collecting higher level taxonomic rank information of closest species match", "\033[0m")))
    #message(cat(paste0("\033[97;", 40, "m","\n", "\033[0m")))
    
    #:::::::::::::::::::::::::::::::::::::
    # Fetch NCBI metadata for sequences
    #:::::::::::::::::::::::::::::::::::::
    id.list <- 1:length(merged.df$NCBI_acc)
    id.chunk <- 250
    fetch.ids <- split(id.list, ceiling(seq_along(id.list)/id.chunk))
    
    fetch.list <- list()
    for (i in 1:length(fetch.ids)) {
      ii <- fetch.ids[[i]]
      while(TRUE){
        r_fetch.try <- try(rentrez::entrez_fetch(db="nucleotide", id=merged.df$NCBI_acc[ii], rettype="gbc")) # OR use rettype="xml" OR useing rettype="gbc", silent=TRUE)
        if(is(r_fetch.try, 'try-error'))
          cat("Failed, trying again in 10 seconds...\n")
        Sys.sleep(3)
        if(!is(r_fetch.try, 'try-error')) break
      }
      r_fetch <- r_fetch.try
      r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_update-date", "INSDSeq_update_date")
      r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_create-date", "INSDSeq_create_date")
      r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_primary.accession", "INSDSeq_primary_accession")
      r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_accession.version", "INSDSeq_accession_version")
      r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_other.seqids", "INSDSeq_other_seqids")
      r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_feature-table", "INSDSeq_feature_table")
      r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_feature-table", "INSDSeq_feature_table")
      #Convert fetched XML file to dataframe format. In this case, the 'INSDSeq' node defines individual records within the XML tree
      fetch.list[[paste("acc", i, sep="_")]] <- lapply(i, function(x) xmlconvert::xml_to_df(text = r_fetch, records.tags = c("//INSDSeq"), check.datatypes=FALSE, check.names=TRUE))
      #svMisc::progress(i, length(fetch.ids))
    }
    
    #Add columns of interest to lookup table
    suppressWarnings(fetch.list.df <- dplyr::bind_rows(fetch.list, .id = "column_label") %>%
                       rowwise() %>%
                       mutate(NCBI_txid = unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'INSDQualifier_name~db_xref||INSDQualifier_value~taxon:', fixed=TRUE))[2], '|', fixed=TRUE))[1]) %>%
                       mutate(isolation_source = unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'isolation_source||INSDQualifier_value~', fixed=TRUE))[2], '|', fixed=TRUE))[1]) %>%
                       mutate(strain = unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'strain||INSDQualifier_value~', fixed=TRUE))[2], '|', fixed=TRUE))[1]) %>%
                       mutate(culture_collection = gsub(":", " ", unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'culture_collection||INSDQualifier_value~', fixed=TRUE))[2], '|', fixed=TRUE))[1])) %>%
                       mutate(culture_collection = ifelse(is.na(culture_collection), strain, culture_collection)) %>%
                       ungroup() %>%
                       mutate(closest_match = paste(stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,1], " ", stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,2], " (", culture_collection, ")", sep="")) %>%
                       mutate(species = paste(stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,1],stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,2], sep=" ")) %>%
                       mutate(taxonomy = gsub(" ", "", .$INSDSeq_taxonomy)) %>%
                       mutate(taxonomy = gsub(";*[Gg]roup*;", ";", taxonomy, perl=TRUE)) %>% #Remove clade rank
                       mutate(taxonomy = gsub(";[A-Za-z]+deae;", ";", taxonomy, perl=TRUE)) %>% #Remove subclass rank
                       mutate(taxonomy = gsub("ales;[A-Za-z]+neae;", "ales;", taxonomy, perl=TRUE)) %>% #Remove suborder ranks
                       mutate(taxonomy = gsub("aceae;[A-Za-z]+eae;", "aceae;", taxonomy, perl=TRUE)) %>% #Remove tribe ranks
                       mutate(rank_domain = stringr::str_split_fixed(.$taxonomy, ";", 7)[,1]) %>%
                       mutate(rank_phylum = stringr::str_split_fixed(.$taxonomy, ";", 7)[,3]) %>%
                       mutate(rank_class = stringr::str_split_fixed(.$taxonomy, ";", 7)[,4]) %>%
                       mutate(rank_order = stringr::str_split_fixed(.$taxonomy, ";", 7)[,5]) %>%
                       mutate(rank_family = stringr::str_split_fixed(.$taxonomy, ";", 7)[,6]) %>%
                       mutate(rank_genus = stringr::str_split_fixed(.$taxonomy, ";", 7)[,7]) %>%
                       mutate(rank_species = gsub("'", "", species)) %>% #Replace instances where single quotations are in species name
                       mutate(genus_tmp = stringr::str_split_fixed(rank_species, " ", 2)[,1]) %>% #Fixing instances where genus is in wrong spot
                       mutate(genus_tmp = gsub("\\[|\\]", "", genus_tmp)) %>% #Fixing instances where genus is in wrong spot
                       mutate(rank_phylum = ifelse(genus_tmp == rank_phylum, "", rank_phylum)) %>% #Fixing instances where genus is in wrong spot
                       mutate(rank_class = ifelse( genus_tmp == rank_class, "", rank_class)) %>% #Fixing instances where genus is in wrong spot
                       mutate(rank_order = ifelse(genus_tmp == rank_order, "", rank_order)) %>% #Fixing instances where genus is in wrong spot
                       mutate(rank_family = ifelse(genus_tmp == rank_family, "", rank_family)) %>% #Fixing instances where genus is in wrong spot
                       mutate(rank_family = ifelse(grepl("aceae$", rank_order), rank_order, rank_family)) %>% #Fixing family upward
                       mutate(rank_family = ifelse(grepl("aceae$", rank_class), rank_class, rank_family)) %>% #Fixing family upward
                       mutate(rank_order = ifelse(rank_order==rank_family, "", rank_order)) %>% #Fixing family upward
                       mutate(rank_class = ifelse(rank_class==rank_family, "", rank_class)) %>% #Fixing family upward
                       mutate(rank_order = ifelse(grepl("ales$", rank_class), rank_class, rank_order)) %>% #Fixing order upward
                       mutate(rank_class = ifelse(rank_class==rank_order, "", rank_class)) %>% #Fixing order upward
                       mutate(species = gsub(" ", "_", species)) %>% #Replace spaces in species name
                       mutate(rank_genus = stringr::str_split_fixed(species, "_", 2)[,1]) %>%
                       mutate(INSDSeq_taxonomy = paste(rank_domain, rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, sep=";")) %>%
                       mutate_at(vars(rank_class, rank_order, rank_family), funs(ifelse(. == "", "NA", .))) %>% #Replace unknown ranks with "NA"
                       mutate(species = paste(stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,1],stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,2], sep=" "))
    )
    
    #::::::::::::::::::::::::::::::::
    if(db=="cpn60"){fetch.list.df$INSDSeq_accession_version <- stringr::str_split_fixed(fetch.list.df$INSDSeq_accession_version, "[.]", 2)[,1]}
    #::::::::::::::::::::::::::::::::
    lookup.list <- setNames(fetch.list.df$INSDSeq_taxonomy, fetch.list.df$INSDSeq_accession_version)
    lookup.list2 <- setNames(fetch.list.df$closest_match, fetch.list.df$INSDSeq_accession_version)
    lookup.list3 <- setNames(fetch.list.df$species, fetch.list.df$INSDSeq_accession_version)
    
    merged.df2 <- merged.df %>%
      mutate(closest_match = dplyr::recode(NCBI_acc, !!!lookup.list2)) %>%
      mutate(taxonomy = dplyr::recode(NCBI_acc, !!!lookup.list)) %>%
      mutate(taxonomy = gsub(" ", "", .$taxonomy)) %>%
      mutate(rank_domain = stringr::str_split_fixed(.$taxonomy, ";", 7)[,1]) %>%
      mutate(rank_phylum = stringr::str_split_fixed(.$taxonomy, ";", 7)[,2]) %>%
      mutate(rank_class = stringr::str_split_fixed(.$taxonomy, ";", 7)[,3]) %>%
      mutate(rank_order = stringr::str_split_fixed(.$taxonomy, ";", 7)[,4]) %>%
      mutate(rank_family = stringr::str_split_fixed(.$taxonomy, ";", 7)[,5]) %>%
      mutate(rank_genus = stringr::str_split_fixed(.$taxonomy, ";", 7)[,6]) %>%
      mutate(rank_species = dplyr::recode(NCBI_acc, !!!lookup.list3)) %>%
      mutate(warning = ifelse(decision =="Fail" | decision =="Pass" & ID <90 & phred_trim <40, "Poor alignment", ""))
  }
  
  if(db=="custom"){
    suppressWarnings(merged.df2 <- merged.df %>%
                       mutate(taxonomy = stringr::str_split_fixed(closest_match, ";", 2)[,2]) %>%
                       mutate(closest_match = stringr::str_split_fixed(closest_match, ";", 8)[,8]) %>% mutate(closest_match = gsub("_", " ", closest_match)) %>%
                       mutate(rank_domain = stringr::str_split_fixed(.$taxonomy, ";", 7)[,1]) %>%
                       mutate(rank_phylum = stringr::str_split_fixed(.$taxonomy, ";", 7)[,2]) %>%
                       mutate(rank_class = stringr::str_split_fixed(.$taxonomy, ";", 7)[,3]) %>%
                       mutate(rank_order = stringr::str_split_fixed(.$taxonomy, ";", 7)[,4]) %>%
                       mutate(rank_family = stringr::str_split_fixed(.$taxonomy, ";", 7)[,5]) %>%
                       mutate(rank_genus = stringr::str_split_fixed(.$taxonomy, ";", 7)[,6]) %>%
                       mutate(rank_species = closest_match) %>%
                       mutate_at(vars(rank_domain, rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, closest_match), funs(gsub("^d__|^p__|^c__|^o__|^f__|^g__|^s  ", "", .))) %>%
                       mutate(warning = ifelse(decision =="Fail" | decision =="Pass" & ID <90 & phred_trim <40, "Poor alignment", ""))
    )
  }
  
  message(cat(paste0("\033[0;", 32, "m","Done.", "\033[0m")))
  
  #::::::::::::::::::::::::
  #Make reactable output
  #::::::::::::::::::::::::
  merged.df3 <- merged.df2 %>% select(warning, date,filename,seqs_raw, phred_raw,Ns_raw, length_raw, phred_spark_raw, seqs_trim, phred_trim,Ns_trim,length_trim, closest_match, NCBI_acc, ID, rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species) %>%
    mutate(rank_genus= paste(stringr::str_split_fixed(rank_genus, ";", 3)[,1])) %>%
    mutate(rank_species = ifelse(grepl("No_match",rank_species) | rank_species == "", "No_match", rank_species)) %>%
    mutate(rank_phylum=ifelse(rank_phylum=="","No_match",rank_phylum)) %>%
    mutate(rank_class=ifelse(rank_class=="","No_match",rank_class)) %>%
    mutate(rank_order=ifelse(rank_order=="","No_match",rank_order)) %>%
    mutate(rank_family=ifelse(rank_family=="","No_match",rank_family)) %>%
    mutate(rank_genus=ifelse(rank_genus=="","No_match",rank_genus)) 
  
  seq.warnings2 <- c()
  seq.warnings2 <- (merged.df3 %>% filter(!warning ==""))$filename
  seq.warnings.txt2 <- paste(seq.warnings2, collapse="\r\n")
  
  if(!length(seq.warnings2)==0){
    message(cat(paste0("\033[97;", 40, "m","\r\nThe following sequences had poor quality alignments and should be manually inspected:","\033[0m","\n     ",
                       "\n", "\033[0;", 95, "m", seq.warnings.txt2,"\033[0m","\n")))
  }
  
  
  #:::::::::::::::::::::::::::::::
  # Outputs
  #:::::::::::::::::::::::::::::::
  
  #Set S4 class for outputs
  isoTAX <- df_to_isoTAX(merged.df3)
  isoTAX@phylum_threshold <- phylum_threshold
  isoTAX@class_threshold <- class_threshold
  isoTAX@order_threshold <- order_threshold
  isoTAX@family_threshold <- family_threshold
  isoTAX@genus_threshold <- genus_threshold
  isoTAX@species_threshold <- species_threshold
  
  #Set dataframe for outputs
  out_df <- S4_to_dataframe(isoTAX)
  
  #####################################################
  
  #Export messages
  message(cat(paste0("\033[97;", 40, "m","'isoTAX' steps completed. Exporting files...", "\033[0m")))
  
  #Export HTML file----------------------------------------------------------------------
  if(export_html == TRUE) {
    output_name <- "02_isoTAX_results.html"
    export_html(isoTAX, quick_search=quick_search, db=db)
    message(cat(paste0("\033[97;", 40, "m","HTML results exported: ", "\033[0m",
                       "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m")))
  }
  #Export BLAST TABLE file----------------------------------------------------------------------
  if(export_blast_table == TRUE) {
    output_name <- "02_isoTAX_results_BLAST_formatted_table.tsv"
    readr::write_tsv(blast.table, file=output_name, col_names=FALSE)
    message(cat(paste0("\033[97;", 40, "m","BLAST format results exported: ", "\033[0m",
                       "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m")))
  }
  #export CSV files----------------------------------------------------------------------
  if(export_csv==TRUE) {
    output_name <- "02_isoTAX_results.csv"
    write.csv(out_df, file=paste(output_name, sep=""), row.names = FALSE)
    message(cat(paste0("\033[97;", 40, "m","CSV file exported:", "\033[0m",
                       "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m",
                       "\033[0;", 31, "m", "  <--- Input for Step 3: 'isoLIB'","\033[0m")))
  }
  
  return(isoTAX)
  ######################################
  
}

