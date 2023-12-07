#' @name isoTAX
#' @title Classify taxonomy of sequences after quality trimming steps.
#' @rdname isoTAX
#' @description  This function performs taxonomic classification steps by searching query Sanger sequences against specified database of interest. Takes CSV input files, extracts FASTA-formatted query sequences and performs global alignment against specified database of interest via Needleman-Wunsch algorithm by wrapping the --usearch_global command implemented in VSEARCH. Default taxonomic rank cutoffs for 16S rRNA gene sequences are based on Yarza et al. 2014, Nat Rev Microbiol.
#' @export
#' @param input Path of CSV output file from isoQC step.
#' @param export_html (Default=TRUE) Output the results as an HTML file
#' @param export_csv (Default=TRUE) Output the results as a CSV file.
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
#' @seealso \code{\link{isoQC}}, \code{\link{isoLIB}}, \code{\link{search_db}}
#' @return Returns taxonomic classification table of class isoTAX. Default taxonomic cutoffs for phylum (75.0), class (78.5), order (82.0), family (86.5), genus (96.5), and species (98.7) demarcation are based on Yarza et al. 2014, Nature Reviews Microbiology (DOI:10.1038/nrmicro3330)
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
                   quick_search=TRUE,
                   db="16S",
                   iddef=2,
                   phylum_cutoff=75.0,
                   class_cutoff=78.5,
                   order_cutoff=82.0,
                   family_cutoff=86.5,
                   genus_cutoff=96.5,
                   species_cutoff=98.7
                   ){


  #Setting file paths-------------------------------------------------
  input <- str_replace_all(input, '\\\\', '/')


  path <- paste(unlist(strsplit(input, '/'))[1:(length(unlist(strsplit(input, '/')))-1)], collapse="/")

  setwd(path)
  vsearch.path.dl <- file.path(system.file("", package="isolateR"), "vsearch")
  suppressWarnings(dir.create(vsearch.path.dl))

    #:::::::::::::::::::::::::::
    #Download VSEARCH software
    #:::::::::::::::::::::::::::
    message(cat(paste0("\n", "\033[97;", 40, "m","Detecting operating system...", "\033[0m")))

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

    #:::::::::::::::::::::::::::::::::::::::::
    #Filtering input file
    #:::::::::::::::::::::::::::::::::::::::::

    isoTAX.input <- read.csv(input)

    #:::::::::::::::::::::::::::::::::::::::::
    #Stage paths for taxonomic classification
    #:::::::::::::::::::::::::::::::::::::::::
    message(cat(paste0("\033[0;", 32, "m","Staging files", "\033[0m")))

    #Stage query FASTA file
    query.fasta.path <- "02_isoTAX_query.fasta"
    make_fasta(csv_file=input,
               col_names= "filename",
               col_seqs= "seqs_trim",
               output=query.fasta.path)

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
                            iddef=iddef)


  message(cat(paste0("\n", "\033[97;", 40, "m","Organizing results", "\033[0m", "\n")))

  #-------------------------------------------------

  query.info <- read.csv(input)
  query.info.list1 <- setNames(query.info$phred_trim, query.info$filename)
  query.info.list2 <- setNames(query.info$date, query.info$filename)

  message(cat(paste0("\033[97;", 40, "m","Merging results...", "\033[0m")))

  merged.df <- merge(query.info, uc.results, by.x="filename", by.y="V9", all=TRUE) %>%
    mutate(NCBI_acc = str_split_fixed(.$V10, " ", 8)[,1]) %>%
    dplyr::rename("ID" = "V4", "closest_match" = "V10") %>%
    mutate(NCBI_acc = ifelse(length_trim < 20, "No_match", NCBI_acc)) %>%
    mutate(closest_match = ifelse(length_trim < 20, "No_match", closest_match)) %>%
    mutate(ID = ifelse(length_trim < 20, 51, ID))

  message(cat(paste0("\033[97;", 40, "m","Done.", "\033[0m")))

  #-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  #------------------------------------------------- library(rentrez)
  message(cat(paste0("\n", "\033[97;", 40, "m","Collecting higher level taxonomic rank information of closest species match", "\033[0m", "\n")))
  message(cat(paste0("\033[97;", 40, "m","\n", "\033[0m")))

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
    svMisc::progress(i, length(fetch.ids))
  }

  message(cat(paste0("\n", "\033[97;", 40, "m","Done.", "\033[0m")))
  message(cat(paste0("\033[97;", 40, "m","Exporting csv file", "\033[0m")))

  #Add columns of interest to lookup table
  fetch.list.df <- as.data.frame(fetch.list$acc_1, stringsasfactors=FALSE) %>%
    rowwise() %>%
    mutate(NCBI_txid = unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'INSDQualifier_name~db_xref||INSDQualifier_value~taxon:', fixed=TRUE))[2], '|', fixed=TRUE))[1]) %>%
    mutate(isolation_source = unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'isolation_source||INSDQualifier_value~', fixed=TRUE))[2], '|', fixed=TRUE))[1]) %>%
    mutate(strain = unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'strain||INSDQualifier_value~', fixed=TRUE))[2], '|', fixed=TRUE))[1]) %>%
    mutate(culture_collection = gsub(":", " ", unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'culture_collection||INSDQualifier_value~', fixed=TRUE))[2], '|', fixed=TRUE))[1])) %>%
    mutate(culture_collection = ifelse(is.na(culture_collection), strain, culture_collection)) %>%
    ungroup() %>%
    mutate(closest_match = paste(str_split_fixed(INSDSeq_organism, " ", 3)[,1], " ", str_split_fixed(INSDSeq_organism, " ", 3)[,2], " (", culture_collection, ")", sep="")) %>%
    mutate(species = paste(str_split_fixed(INSDSeq_organism, " ", 3)[,1],str_split_fixed(INSDSeq_organism, " ", 3)[,2], sep=" "))

  #::::::::::::::::::::::::::::::::
  lookup.list <- setNames(fetch.list.df$INSDSeq_taxonomy, fetch.list.df$INSDSeq_accession_version)
  lookup.list2 <- setNames(fetch.list.df$closest_match, fetch.list.df$INSDSeq_accession_version)
  lookup.list3 <- setNames(fetch.list.df$species, fetch.list.df$INSDSeq_accession_version)


  merged.df2 <- merged.df %>%
    mutate(closest_match = dplyr::recode(NCBI_acc, !!!lookup.list2)) %>%
    mutate(taxonomy = dplyr::recode(NCBI_acc, !!!lookup.list)) %>%
    mutate(taxonomy = gsub(" ", "", .$taxonomy)) %>%
    mutate(rank_domain = str_split_fixed(.$taxonomy, ";", 6)[,1]) %>%
    mutate(rank_phylum = str_split_fixed(.$taxonomy, ";", 6)[,2]) %>%
    mutate(rank_class = str_split_fixed(.$taxonomy, ";", 6)[,3]) %>%
    mutate(rank_order = str_split_fixed(.$taxonomy, ";", 6)[,4]) %>%
    mutate(rank_family = str_split_fixed(.$taxonomy, ";", 6)[,5]) %>%
    mutate(rank_genus = str_split_fixed(.$taxonomy, ";", 6)[,6]) %>%
    mutate(rank_species = dplyr::recode(NCBI_acc, !!!lookup.list3)) %>%
    mutate(warning = ifelse(decision =="Fail" | decision =="Pass" & ID <90 & phred_trim <40, "Poor alignment", ""))


  message(cat(paste0("\n", "\033[97;", 40, "m","Done.", "\033[0m")))

  #::::::::::::::::::::::::
  #Make reactable output
  #::::::::::::::::::::::::
  merged.df3 <- merged.df2 %>% select(warning, date,filename,seqs_raw, phred_raw,Ns_raw, length_raw, phred_spark_raw, seqs_trim, phred_trim,Ns_trim,length_trim, closest_match, NCBI_acc, ID, rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species) %>%
    mutate(rank_genus= paste(str_split_fixed(rank_genus, ";", 3)[,1])) %>%
    mutate(rank_species = ifelse(grepl("No_match",rank_species) | rank_species == "", "No_match", rank_species)) %>%
    mutate(rank_phylum=ifelse(rank_phylum=="","No_match",rank_phylum)) %>%
    mutate(rank_class=ifelse(rank_class=="","No_match",rank_class)) %>%
    mutate(rank_order=ifelse(rank_order=="","No_match",rank_order)) %>%
    mutate(rank_family=ifelse(rank_family=="","No_match",rank_family)) %>%
    mutate(rank_genus=ifelse(rank_genus=="","No_match",rank_genus)) %>%
    mutate(filename = str_remove(filename,'_[:alpha:][:digit:][:digit:]\\.ab1$'))

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
  isoTAX@phylum_cutoff <- phylum_cutoff
  isoTAX@class_cutoff <- class_cutoff
  isoTAX@order_cutoff <- order_cutoff
  isoTAX@family_cutoff <- family_cutoff
  isoTAX@genus_cutoff <- genus_cutoff
  isoTAX@species_cutoff <- species_cutoff

  #Set dataframe for outputs
  out_df <- S4_to_dataframe(isoTAX)

  #####################################################

  #Export messages
  message(cat(paste0("\033[97;", 40, "m","'isoTAX' steps completed. Exporting files...", "\033[0m")))

  #Export HTML file----------------------------------------------------------------------
  if(export_html == TRUE) {
    output_name <- "02_isoTAX_results.html"
    export_html(isoTAX)
    #merged.df3_react <- export_html2(out_df, path, output=output_name)
    message(cat(paste0("\033[97;", 40, "m","HTML results exported: ", "\033[0m",
                       "\033[0;", 32, "m", " ", file.path(path, output_name),"\033[0m", "\n")))
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

