

assign_taxonomy <- function(folder=NULL,
                            export_csv=TRUE,
                            verbose=TRUE,
                            exclude=NULL,
                            skip_search=FALSE,
                            quick_search=FALSE,
                            add_fungi_db=FALSE){
  
  # function requirements------------------------------------------------------------
#checking for required packages; installing those not yet installed
  suppressWarnings({
    if(require(dplyr)==FALSE) install.packages('dplyr')
if(require(stringr)==FALSE) install.packages('stringr')
if(require(R.utils)==FALSE) install.packages('R.utils')
if(require(rentrez)==FALSE) install.packages('rentrez')
if(require(xmlconvert)==FALSE) install.packages('xmlconvert')
if(require(reactable)==FALSE) install.packages('reactable')
if(require(reactablefmtr)==FALSE) install.packages('reactablefmtr')
if(require(pander)==FALSE) install.packages('pander')
if(require(svMisc)==FALSE) install.packages('svMisc')
if(require(Biostrings)==FALSE){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("Biostrings", update = FALSE)
  }
if(require(seqinr)==FALSE) install.packages('seqinr')
  
#loading required packages
library(dplyr)
library(stringr)
library(R.utils)
library(rentrez)
library(xmlconvert)
library(reactable)
library(reactablefmtr)
library(pander)
library(svMisc)
library(Biostrings)
library(seqinr)
  })
  
  
#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------
#-------------------------------------------------

#Setting file paths-------------------------------------------------

path <- str_replace_all(folder, '\\\\', '/')
suppressWarnings(dir.create(file.path(path, 'output/NCBI_databases')))
files  <- dir(file.path(path, 'output/NCBI_databases'), full.names = TRUE)

#------------------------------------------------- functions
df_to_fasta <- function(df= NULL, name_col = NULL, seq_col=NULL ){
  fasta.out <- Biostrings::DNAStringSet(setNames(df[,paste(seq_col)], df[,paste(name_col)]))
  return(fasta.out)
}
#------------------------------------------------- 

if(skip_search==FALSE){
  
#!file.exists(file.path(path, 'output/NCBI_databases/bacteria.16SrRNA.fna.gz'))
  #!file.exists(file.path(path, 'output/NCBI_databases/bacteria.16SrRNA.fna.gz'))
  
#:::::::::::::::::::::::::::::::::
#Download NCBI 16S rRNA database
#:::::::::::::::::::::::::::::::::
  
  if(add_fungi_db == TRUE){
    download.file("https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz", file.path(path, 'output/NCBI_databases/bacteria.16SrRNA.fna.gz'), mode='wb')
    R.utils::gunzip(file.path(path,"output/NCBI_databases/bacteria.16SrRNA.fna.gz"), remove = FALSE, overwrite=TRUE)
    message(cat(paste0("\n", "\033[0;", 32, "m","Download complete for NCBI Bacterial 16S rRNA database.", "\033[0m", "\n")))
    db.fasta <- "output/NCBI_databases/bacteria.16SrRNA.fna"
    message(cat(paste0("\n", "\033[0;", 32, "m","Adding ITS sequences.", "\033[0m", "\n")))
    download.file("https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz", file.path(path, 'output/NCBI_databases/fungi.ITS.fna.gz'), mode='wb')
    R.utils::gunzip(file.path(path,"output/NCBI_databases/fungi.ITS.fna.gz"), remove = FALSE, overwrite=TRUE)
    message(cat(paste0("\n", "\033[0;", 32, "m","Download complete for NCBI Fungal ITS database.", "\033[0m", "\n")))
    #Concatenate 16S and ITS dbs
    db.concat <- c(seqinr::read.fasta(file=file.path(path,"output/NCBI_databases/bacteria.16SrRNA.fna"), as.string=TRUE, forceDNAtolower=FALSE),
                   seqinr::read.fasta(file=file.path(path,"output/NCBI_databases/fungi.ITS.fna"), as.string=TRUE, forceDNAtolower=FALSE))
    #Write concatenated db to fasta
    seqinr::write.fasta(sequences=db.concat,
                        names=names(db.concat),
                        as.string=TRUE, nbchar = 1000, file.out=file.path(path,"output/NCBI_databases/16S_ITS_concat.fna"))
    #Set path
    db.fasta <- file.path("output/NCBI_databases/16S_ITS_concat.fna")
  }  else {
    download.file("https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz", file.path(path, 'output/NCBI_databases/bacteria.16SrRNA.fna.gz'), mode='wb')
    R.utils::gunzip(file.path(path,"output/NCBI_databases/bacteria.16SrRNA.fna.gz"), remove = FALSE, overwrite=TRUE)
    message(cat(paste0("\n", "\033[0;", 32, "m","Download complete for NCBI Bacterial 16S rRNA database.", "\033[0m", "\n")))
    db.fasta <- file.path("output/NCBI_databases/bacteria.16SrRNA.fna")
  }
  

gz_files <- stringr::str_subset(files, 'vsearch')


#:::::::::::::::::::::::::::
#Download USEARCH software
#:::::::::::::::::::::::::::
message(cat(paste0("\n", "\033[97;", 40, "m","Detecting operating system...", "\033[0m", "\n")))

vsearch_files <- stringr::str_subset(files, 'vsearch')

get_os <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "osx-mac"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "osx-mac"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

if(paste(get_os())=="windows"){
  if(identical(vsearch_files, character(0))){
    message(cat(paste0("\n", "\033[0;", 32, "m","Operating system is ---> Windows-based <---", "\033[0m", "\n")))
    download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-win-x86_64.zip", file.path(path, 'output/NCBI_databases/vsearch-2.23.0-win-x86_64.zip'), mode='wb')
    unzip(file.path(path,"output/NCBI_databases/vsearch-2.23.0-win-x86_64.zip"),  exdir=file.path(path,"output/NCBI_databases"))
    file.copy(file.path(path, "output/NCBI_databases/vsearch-2.23.0-win-x86_64/bin/vsearch.exe"), file.path(path, "output/NCBI_databases/vsearch-2.23.0.exe"), overwrite=TRUE)
    unlink(file.path(path,"output/NCBI_databases/vsearch-2.23.0-win-x86_64"),recursive=TRUE)
    message(cat(paste0("\n", "\033[0;", 32, "m","Download complete.", "\033[0m", "\n")))
    vsearch.path <- file.path(path,"output/NCBI_databases/vsearch-2.23.0.exe")
  } else {
    message(cat(paste0("\n", "\033[0;", 32, "m","Operating system is ---> Windows-based <---", "\033[0m")))
    message(cat(paste0("\033[0;", 32, "m","No additional download needed.", "\033[0m", "\n")))
    vsearch.path <- file.path(path,"output/NCBI_databases/vsearch-2.23.0.exe")
  }
}

if(paste(get_os())=="osx-mac"){
  if(identical(vsearch_files, character(0))){
  message(cat(paste0("\n", "\033[0;", 32, "m","Operating system is ---> MacOS-based <---", "\033[0m", "\n")))
  download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-macos-universal.tar.gz", file.path(path, 'output/NCBI_databases/vsearch-2.23.0-macos-universal.tar.gz'), mode='wb')
  #R.utils::gunzip(file.path(path,"output/NCBI_databases/vsearch-2.23.0-macos-universal.tar.gz"), remove = FALSE, overwrite=TRUE)
  untar(file.path(path,"output/NCBI_databases/vsearch-2.23.0-macos-universal.tar.gz"),  exdir=file.path(path,"output"))
  file.copy(file.path(path, "output/NCBI_databases/vsearch-2.23.0-macos-universal/bin/vsearch"), file.path(path, "output/NCBI_databases/vsearch-2.23.0_macos"), overwrite=TRUE)
  unlink(file.path(path,"output/NCBI_databases/vsearch-2.23.0-macos-universal"),recursive=TRUE)
  message(cat(paste0("\n", "\033[0;", 32, "m","Download complete.", "\033[0m", "\n")))
  vsearch.path <- file.path(path,"output/NCBI_databases/vsearch-2.23.0_macos")
  } else {
    message(cat(paste0("\n", "\033[0;", 32, "m","Operating system is ---> MacOS-based <---", "\033[0m")))
    message(cat(paste0("\033[0;", 32, "m","No additional download needed.", "\033[0m", "\n")))
    vsearch.path <- file.path(path,"output/NCBI_databases/vsearch-2.23.0_macos")
  }
}

if(paste(get_os())=="linux"){
  if(identical(vsearch_files, character(0))){
    message(cat(paste0("\n", "\033[0;", 32, "m","Operating system is ---> Linux-based <---", "\033[0m", "\n")))
    download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-linux-x86_64.tar.gz", file.path(path, 'output/NCBI_databases/vsearch-2.23.0-linux-x86_64.tar.gz'), mode='wb')
    untar(file.path(path,"output/NCBI_databases/vsearch-2.23.0-linux-x86_64.tar.gz"),  exdir=file.path(path,"output"))
    file.copy(file.path(path, "output/NCBI_databases/vsearch-2.23.0-linux-x86_64/bin/vsearch"), file.path(path, "output/NCBI_databases/vsearch-2.23.0"), overwrite=TRUE)
    unlink(file.path(path,"output/NCBI_databases/vsearch-2.23.0-linux-x86_64"),recursive=TRUE)
    message(cat(paste0("\n", "\033[0;", 32, "m","Download complete.", "\033[0m", "\n")))
    vsearch.path <- file.path(path,"output/NCBI_databases/vsearch-2.23.0")
  } else {
    message(cat(paste0("\n", "\033[0;", 32, "m","Operating system is ---> Linux-based <---", "\033[0m")))
    message(cat(paste0("\033[0;", 32, "m","No additional download needed.", "\033[0m", "\n")))
    vsearch.path <- file.path(path,"output/NCBI_databases/vsearch-2.23.0")
  }
  
}

#::::::::::::
#Unzip files
#::::::::::::

#gz_files <- stringr::str_subset(files, '.gz')
#if(! identical(gz_files, character(0))){
#unlink("output/NCBI_databases/*.tmp", recursive = TRUE, force = FALSE)
#message(cat(paste0("\n", "\033[0;", 32, "m","Unzipping files...", "\033[0m", "\n")))
#lapply(gz_files, R.utils::gunzip, remove = FALSE, overwrite=TRUE)
#message(cat(paste0("\n", "\033[0;", 32, "m","Done", "\033[0m", "\n")))
#}

#::::::::::::
#Stage paths
#::::::::::::
message(cat(paste0("\n", "\033[0;", 32, "m","Staging files", "\033[0m", "\n")))

setwd(path)
#query.in <- readDNAStringSet(paste0("abif_fasta2_output___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))],".fasta", sep=""), format="fasta")
#query.df <- as.data.frame(cbind(names(query.in), paste(query.in))) %>% mutate(length= nchar(V2)) %>% mutate(V2 = ifelse(length < 200, "NA", V2))
#write.fasta(sequences=as.list(query.df$V2), names=query.df$V1, as.string=TRUE, nbchar = 1000, file.out="query.fasta")
#query.fasta <- "query.fasta"

query.fasta <- paste0("output/abif_fasta2_output_PASS___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))],".fasta", sep="")
#query.fasta <- "V3-V6seq.FASTA"
#db.fasta <- "output/NCBI_databases/bacteria.16SrRNA.fna"
#db.fasta <- db.fasta
uc.out <- paste0("output/assign_taxonomy_output_PASS___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))],".uc", sep="")
b6.out <- paste0("output/assign_taxonomy_output_PASS___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))],".b6o", sep="")
#db.fasta <- file.path(path, "output/bacteria.16SrRNA.fna")
#uc.out <- file.path(path, "V3-V6seq_results.uc")
usearch_files <- stringr::str_subset(files, 'usearch')
usearch.path <- gsub(".gz", "", usearch_files[1], fixed=TRUE)

#:::::::::::::::::::::
#Search NCBI database
#:::::::::::::::::::::
message(cat(paste0("\n", "\033[97;", 40, "m","Searching query sequences against NCBI database", "\033[0m")))
message(cat(paste0("\033[0;", 32, "m","This step will take 2-3 minutes...", "\033[0m", "\n")))

#system2(usearch.path, paste(" --usearch_global ", query.fasta, " --db ", db.fasta, " --output_no_hits --blast6out ", b6.out, " --uc ", uc.out, " --id 0.7 --maxaccepts 0 --maxrejects 0 --top_hits_only --strand both --threads 1 ", sep=""), stdout="", stderr="")
#system2(vsearch.path, paste(" --usearch_global ", query.fasta, " --db ", db.fasta, " --output_no_hits --blast6out ", b6.out, " --uc ", uc.out, " --id 0.7 --maxaccepts 0 --maxrejects 0 --top_hits_only --strand both --threads 1 ", sep=""), stdout="", stderr="")
if(quick_search==TRUE){
  system2(vsearch.path, paste(" --usearch_global ", query.fasta, " --db ", db.fasta, " --blast6out ", b6.out, " --uc ", uc.out, " --id 0.7 --maxaccepts 100 --maxrejects 100 --top_hits_only --strand both", sep=""), stdout="", stderr="")
}

if(quick_search==FALSE){
  system2(vsearch.path, paste(" --usearch_global ", query.fasta, " --db ", db.fasta, " --blast6out ", b6.out, " --uc ", uc.out, " --id 0.7 --maxaccepts 0 --maxrejects 0 --top_hits_only --strand both ", sep=""), stdout="", stderr="")
}

} #end of "skip_search=TRUE"

#:::::::::::::::::
#Organize results
#:::::::::::::::::
message(cat(paste0("\n", "\033[97;", 40, "m","Determining closest species match", "\033[0m", "\n")))

uc.results <- read.csv(uc.out, sep="\t", header = FALSE) %>% 
  mutate(V4 = gsub("*", "51", .$V4, fixed=TRUE)) %>% mutate_at(vars(V4), as.numeric) %>%   # Fix percent ID column to make numeric
  mutate(V10 = gsub("*", "No_match", .$V10, fixed=TRUE)) %>%                               # Fix unmatched taxa column
  arrange(desc(V4)) %>% distinct(V9, .keep_all=TRUE) #%>% select(V1, V4, V9, V10) #%>% filter(V10 != "*")

message(cat(paste0("\n", "\033[97;", 40, "m","Organizing results", "\033[0m", "\n")))

#-------------------------------------------------
query.seqs <- Biostrings::readBStringSet(query.fasta)
query.seqs.df <- as.data.frame(names(query.seqs)) %>% mutate(V9 = names(query.seqs), query_seq = paste(query.seqs)) %>% mutate(length= nchar(query_seq)) %>% mutate(Ns = str_count(.$query_seq, "[Nn]")) %>% select(-1)
query.info <- read.csv(paste("output/abif_fasta2_output_PASS___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))],".csv", sep=""), row.names= 1)
query.info.list1 <- setNames(query.info$phred_trim, query.info$filename)
query.info.list2 <- setNames(query.info$date, query.info$filename)


message(cat(paste0("\n", "\033[97;", 40, "m","Merging results...", "\033[0m", "\n")))

merged.df <- merge(uc.results, query.seqs.df, by="V9", all=TRUE) %>% arrange(V9) %>%
  mutate(NCBI_acc = str_split_fixed(.$V10, " ", 8)[,1]) %>%
  mutate(species = paste(str_split_fixed(.$V10, " ", 8)[,2], str_split_fixed(.$V10, " ", 8)[,3], sep=" ")) %>%
  dplyr::rename("ID" = "V4", "filename" = "V9", "closest_match" = "V10") %>%
  mutate(phred_trim = dplyr::recode(filename, !!!query.info.list1)) %>%
  mutate(date = dplyr::recode(filename, !!!query.info.list2)) %>%
  #mutate(ID = ifelse(length <200 | phred_trim < 20, 51.0, ID)) %>%
  #mutate(NCBI_acc = ifelse(length <200 | phred_trim < 20, "No_match", NCBI_acc)) %>%
  #mutate(closest_match = ifelse(length <200 | phred_trim < 20, "No_match", NCBI_acc)) %>%
  select(date,filename,species,ID,phred_trim,length,Ns,NCBI_acc,closest_match,query_seq)


message(cat(paste0("\n", "\033[97;", 40, "m","Done.", "\033[0m", "\n")))

#write.csv(merged.df, file=paste("V3-V6seq_results_", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))],".csv", sep=""))

#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------- library(rentrez)
message(cat(paste0("\n", "\033[97;", 40, "m","Collecting higher level taxonomic rank information of closest species match", "\033[0m", "\n", "\r\n")))
message(cat(paste0("\n", "\033[97;", 40, "m","\n", "\033[0m", "\n")))

#::::::::::::::::::::::::::::::::
# Step 1b) For loop
#::::::::::::::::::::::::::::::::
id.list <- 1:length(merged.df$NCBI_acc)
id.chunk <- 250
fetch.ids <- split(id.list, ceiling(seq_along(id.list)/id.chunk))

fetch.list <- list()
for (i in 1:length(fetch.ids)) {
  #setTxtProgressBar(pb.bar,i)
  #t1 <- Sys.time()
  ii <- fetch.ids[[i]]
  while(TRUE){
    r_fetch.try <- try(entrez_fetch(db="nucleotide", id=merged.df$NCBI_acc[ii], rettype="gbc")) # OR use rettype="xml" OR useing rettype="gbc", silent=TRUE)
    if(is(r_fetch.try, 'try-error'))
      cat("Failed, trying again in 10 seconds...\n")
    Sys.sleep(3)
    if(!is(r_fetch.try, 'try-error')) break
  }
  r_fetch <- r_fetch.try
  #r_fetch <- entrez_fetch(db="nucleotide", id=merged.df$NCBI_acc[ii], rettype="gbc") # OR use rettype="xml" OR useing rettype="gbc"
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_update-date", "INSDSeq_update_date")
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_create-date", "INSDSeq_create_date")
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_primary.accession", "INSDSeq_primary_accession")
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_accession.version", "INSDSeq_accession_version")
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_other.seqids", "INSDSeq_other_seqids")
  r_fetch <- stringr::str_replace_all(r_fetch, "INSDSeq_feature-table", "INSDSeq_feature_table")
  #Convert fetched XML file to dataframe format. In this case, the 'INSDSeq' node defines individual records within the XML tree
  #xml.df <- xml_to_df(text = r_fetch, records.tags = c("//INSDSeq"), check.datatypes=FALSE, check.names=TRUE)
  fetch.list[[paste("acc", i, sep="_")]] <- lapply(i, function(x) xmlconvert::xml_to_df(text = r_fetch, records.tags = c("//INSDSeq"), check.datatypes=FALSE, check.names=TRUE))
  #fetch.list[[paste("acc_1")]] <- rbind(as.data.frame(fetch.list[[1]]), as.data.frame(fetch.list[[paste("acc", i, sep="_")]]))
  #fetch.list[[paste("acc", i, sep="_")]] <- NULL
  svMisc::progress(i, length(fetch.ids))
  #svMisc::progress(((i/length(fetch.ids))*100), 100, progress.bar=TRUE)
  #Sys.sleep(0.01)
  #t2 <- Sys.time()
  #print(difftime(t2, t1, units = "secs")[[1]])
  #if (((i/length(fetch.ids))*100) == 100) message("Done!")
}
message(cat(paste0("\n", "\033[97;", 40, "m","Done.", "\033[0m", "\n")))

message(cat(paste0("\n", "\033[97;", 40, "m","Exporting csv file", "\033[0m", "\n")))

#fetch.list.df <- bind_rows(fetch.list, .id = "column_label") #Only works if more than 2 entries in the list
fetch.list.df <- as.data.frame(fetch.list$acc_1, stringsasfactors=FALSE) #Use this if only 1 entry in the list

#::::::::::::::::::::::::::::::::
lookup.list <- setNames(fetch.list.df$INSDSeq_taxonomy, fetch.list.df$INSDSeq_accession_version)
lookup.list2 <- setNames(fetch.list.df$INSDSeq_organism, fetch.list.df$INSDSeq_accession_version)

merged.df2 <- merged.df %>% mutate(taxonomy = dplyr::recode(NCBI_acc, !!!lookup.list)) %>%
  mutate(species = dplyr::recode(NCBI_acc, !!!lookup.list2)) %>%
  mutate(taxonomy = gsub(" ", "", .$taxonomy)) %>%
  mutate(domain = str_split_fixed(.$taxonomy, ";", 6)[,1]) %>%
  mutate(phylum = str_split_fixed(.$taxonomy, ";", 6)[,2]) %>%
  mutate(class = str_split_fixed(.$taxonomy, ";", 6)[,3]) %>%
  mutate(order = str_split_fixed(.$taxonomy, ";", 6)[,4]) %>%
  mutate(family = str_split_fixed(.$taxonomy, ";", 6)[,5]) %>%
  mutate(genus = str_split_fixed(.$taxonomy, ";", 6)[,6]) %>%
  mutate(warning = ifelse(ID <90 & phred_trim <45, "Poor alignment", ""))


message(cat(paste0("\r\n\r\n", "\033[97;", 40, "m","Done.", "\033[0m", "\n")))

#::::::::::::::::::::::::
#Make reactable output
#::::::::::::::::::::::::


merged.df3 <- merged.df2 %>% select(warning,date,filename,ID,species,NCBI_acc, phred_trim, Ns,length,query_seq,phylum,class,order,family,genus) %>%
                             mutate(species2 = species) %>%
                             mutate(genus= paste(str_split_fixed(genus, ";", 3)[,1])) %>%
                             mutate(species2= paste(str_split_fixed(species2, " ", 3)[,1],str_split_fixed(species2, " ", 3)[,2], sep=" ")) %>%
                             mutate(species2 = ifelse(grepl("No_match",species2) | species2 == "", "No_match", species2)) %>%
                             mutate(species = species2) %>%
                             #mutate(genus = str_split_fixed(genus, ";", 2)[,1]) %>%
                             #mutate(species2 = gsub("_|No_match", "", .$species2, fixed=TRUE)) %>% 
                             mutate(seq_warnings_col = ifelse(ID <90 & phred_trim <45, grDevices::adjustcolor("#B20A2C", alpha.f=0.7), NA)) %>%
                             mutate(species_colors = dplyr::case_when(genus != "" ~ "green", TRUE ~ "red")) %>% 
                             mutate(species_col = ifelse(ID > 98.9, "#31a354", NA)) %>%
                             mutate(genus_col = ifelse(ID > 95.7, "#31a354", NA)) %>%
                             mutate(family_col = ifelse(ID > 89.5, "#31a354", NA)) %>%
                             mutate(order_col = ifelse(ID > 86.1, "#31a354", NA)) %>%
                             mutate(class_col = ifelse(ID > 83.4, "#31a354", NA)) %>%
                             mutate(phylum_col = ifelse(ID > 80.7, "#31a354", NA)) %>%
                             mutate(phylum=ifelse(phylum=="","No_match",phylum)) %>%
                             mutate(class=ifelse(class=="","No_match",class)) %>%
                             mutate(order=ifelse(order=="","No_match",order)) %>%
                             mutate(family=ifelse(family=="","No_match",family)) %>%
                             mutate(genus=ifelse(genus=="","No_match",genus))
  

seq.warnings2 <- c()
seq.warnings2 <- (merged.df3 %>% filter(ID < 90 & phred_trim <45))$filename
seq.warnings.txt2 <- paste(seq.warnings2, collapse="\r\n     ")

message(cat(paste0("\033[97;", 40, "m","\r\nThe following sequences had poor quality alignments and should be manually inspected before next steps:","\033[0m","\n     ",
                   "\033[0;", 95, "m", seq.warnings.txt2,"\033[0m","\n")))

#merged.df3 <- merged.df3 %>% filter(!filename %in% seq.warnings2)

#grDevices::adjustcolor("#31a354", alpha.f = 0.5) #97D1A9
#grDevices::adjustcolor(merged.df3$genus_col[[index]], alpha.f = 0.5)
#which(!is.na(merged.df3$species_col))
#merged.df3$species_colors

#:::::::::::::::::::::::::::::::
#List function for table below
#:::::::::::::::::::::::::::::::
list_with_names<-function(...){dplyr::lst(...)}
uni.font <- 9


#::::::::::::::::::::::::
#Reactable
#::::::::::::::::::::::::
merged.df3_react <- reactable(merged.df3,
                            fullWidth=FALSE, searchable = TRUE, bordered = TRUE, resizable =TRUE, #width = 1200,
                            defaultPageSize=200, highlight = TRUE, showSortable = TRUE, compact=TRUE, wrap = TRUE, #highlight=TRUE,
                            #theme = reactableTheme(headerStyle = list(wrap=TRUE)),
                            #theme = fivethirtyeight(centered = TRUE, header_font_size = 11),
                            theme = reactableTheme(                  #borderColor = "#555",
                              headerStyle = list(
                                fontFamily = "sans-serif", fontWeight="bold", fontSize=uni.font,
                               "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                               "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"))),
                            #borderColor = "#555")),
                            #defaultColDef = colDef(minWidth = 100), #
                            #defaultColDef = colDef(style = reactablefmtr::cell_style(checkseq.sub, font_size=12)),
                            style = list(minWidth = 1400),
                            #--------------------Universal settings
                            defaultColDef = colDef(vAlign = "center", style = list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                            #--------------------Column specific settings
                            columns = list(
                                species2 = colDef(name="species", vAlign = "center", headerVAlign = "bottom", minWidth=130, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                               style = function(value, index, species2){
                                                 if(is.na(merged.df3$species_col[[index]])) { 
                                                   color.x <- "black" 
                                                 } else if(!is.na(merged.df3$species_col[[index]])) {
                                                   color.x <- "white" }
                                                 list(background = merged.df3$species_col[[index]], color = color.x, 
                                                      whiteSpace = "nowrap",
                                                      fontSize=uni.font,
                                                      #fontWeight="bold",
                                                      align="center",
                                                      fontFamily = "sans-serif")}),
                                genus = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=100, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                style = function(value, index, genus){
                                                  if(is.na(merged.df3$genus_col[[index]])) { 
                                                    color.x <- "black" 
                                                  } else if(!is.na(merged.df3$genus_col[[index]])) {
                                                    color.x <- "white" }
                                                  list(background = merged.df3$genus_col[[index]], color = color.x, 
                                                       whiteSpace = "nowrap",
                                                       fontSize=uni.font,
                                                       #fontWeight="bold",
                                                       align="center",
                                                       fontFamily = "sans-serif")}),
                                family = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=100, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                               style = function(value, index, family){
                                                 if(is.na(merged.df3$family_col[[index]])) { 
                                                   color.x <- "black" 
                                                 } else if(!is.na(merged.df3$family_col[[index]])) {
                                                   color.x <- "white" }
                                                  list(background = merged.df3$family_col[[index]], color = color.x, 
                                                       whiteSpace = "nowrap",
                                                       fontSize=uni.font, 
                                                       #fontWeight="bold",
                                                       align="center",
                                                       fontFamily = "sans-serif")}),
                                order = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=100, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                               style = function(value, index, order){
                                                 if(is.na(merged.df3$order_col[[index]])) { 
                                                   color.x <- "black" 
                                                 } else if(!is.na(merged.df3$order_col[[index]])) {
                                                   color.x <- "white" }
                                                 list(background = merged.df3$order_col[[index]], color = color.x, 
                                                      whiteSpace = "nowrap",
                                                      fontSize=uni.font, 
                                                      #fontWeight="bold",
                                                      align="center",
                                                      fontFamily = "sans-serif")}),
                                class = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=100, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                style = function(value, index, class){
                                                  if(is.na(merged.df3$class_col[[index]])) { 
                                                    color.x <- "black" 
                                                  } else if(!is.na(merged.df3$class_col[[index]])) {
                                                    color.x <- "white" }
                                                  list(background = merged.df3$class_col[[index]], color = color.x, 
                                                       whiteSpace = "nowrap",
                                                       fontSize=uni.font, 
                                                       #fontWeight="bold",
                                                       align="center",
                                                       fontFamily = "sans-serif")}),
                                phylum = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=100, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                               style = function(value, index, phylum){
                                                 if(is.na(merged.df3$phylum_col[[index]])) { 
                                                   color.x <- "black" 
                                                 } else if(!is.na(merged.df3$phylum_col[[index]])) {
                                                   color.x <- "white" }
                                                 list(background = merged.df3$phylum_col[[index]], color = color.x, 
                                                      whiteSpace = "nowrap",
                                                      fontSize=uni.font, 
                                                      #fontWeight="bold",
                                                      align="center",
                                                      fontFamily = "sans-serif")}),
                                warning = colDef(name="Warning", vAlign = "center", headerVAlign = "bottom", width=75, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                 style = function(value, index, alignment){
                                                   if(is.na(merged.df3$seq_warnings_col[[index]])) { 
                                                     color.x <- "black" 
                                                   } else if(!is.na(merged.df3$seq_warnings_col[[index]])) {
                                                     color.x <- "white" }
                                                   list(background = merged.df3$seq_warnings_col[[index]], color = color.x, 
                                                        whiteSpace = "nowrap",
                                                        fontSize=uni.font, 
                                                        #fontWeight="bold",
                                                        align="center",
                                                        fontFamily = "sans-serif")}),
                                filename = colDef(vAlign = "center", headerVAlign = "bottom", minWidth=150, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                  style = function(value, index, filename){
                                                    if(is.na(merged.df3$seq_warnings_col[[index]])) { 
                                                      color.x <- "black" 
                                                    } else if(!is.na(merged.df3$seq_warnings_col[[index]])) {
                                                      color.x <- "white" }
                                                    list(background = merged.df3$seq_warnings_col[[index]], color = color.x, 
                                                         whiteSpace = "nowrap",
                                                         fontSize=uni.font, 
                                                         #fontWeight="bold",
                                                         align="center",
                                                         fontFamily = "sans-serif")}),
                                date = colDef(vAlign = "center", headerVAlign = "bottom", width=65, #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                                  style = function(value, index, date){
                                                    if(is.na(merged.df3$seq_warnings_col[[index]])) { 
                                                      color.x <- "black" 
                                                    } else if(!is.na(merged.df3$seq_warnings_col[[index]])) {
                                                      color.x <- "white" }
                                                    list(background = merged.df3$seq_warnings_col[[index]], color = color.x, 
                                                         whiteSpace = "nowrap",
                                                         fontSize=uni.font, 
                                                         #fontWeight="bold",
                                                         align="center",
                                                         fontFamily = "sans-serif")}),
                                species = colDef(name="Closest match", minWidth=170), #style = list(fontWeight="bold")),# cell = color_tiles(merged.df3, color_ref = 'species_col')),
                                NCBI_acc = colDef(width=70),
                                query_seq = colDef(minWidth=70),
                                seq_warnings_col = colDef(show=FALSE),
                                species_colors = colDef(show=FALSE),
                                species_col = colDef(show=FALSE),
                                genus_col = colDef(show=FALSE),
                                family_col = colDef(show=FALSE),
                                class_col = colDef(show=FALSE),
                                order_col = colDef(show=FALSE),
                                phylum_col = colDef(show=FALSE),
                                # Gauge widgets
                                #-------------------
                                ID = colDef(name = "ID", width= 60, 
                                            cell = gauge_chart(merged.df3, #fill_color_ref = "ID_log",
                                                   #fill_color = '#1A9641',
                                                   #fill_color = c('#FDAE61', '#FFFFBF','#A6D96A','#1A9641'), opacity = 0.5,
                                                   fill_color = c('#e5f5e0', '#a1d99b', '#31a354'), opacity = 0.5,
                                                   #viridis::viridis(5)[1:4]),
                                                   #'#D7191C','#FDAE61','#FFFFBF','#A6D96A','#1A9641'),
                                                   #number_fmt = scales::number_format(scale=50), #scales::log_trans(), #scales::comma,
                                                   #bold_text = TRUE,
                                                   text_size = 10,
                                                   #background = '#555555',
                                                   min_value = 51, max_value=100,
                                                   #max_value = max(c(checkseq.sub$phred_raw,checkseq.sub$phred_trim)),
                                                   show_min_max = FALSE)),
                                phred_trim = colDef(name = "trimQ", width= 60, cell = gauge_chart(merged.df3,
                                                     fill_color = c('#e5f5e0', '#a1d99b', '#31a354'), opacity = 0.5,
                                                     number_fmt = scales::comma,
                                                     text_size = 10,
                                                     max_value = max(merged.df3$phred_trim),
                                                     show_min_max = FALSE)),
                                Ns = colDef(name = "Ns", width= 60, 
                                            cell = gauge_chart(merged.df3, #fill_color_ref = "ID_log",
                                                               #fill_color = '#1A9641',
                                                               #fill_color = c('#FDAE61', '#FFFFBF','#A6D96A','#1A9641'), opacity = 0.5,
                                                               fill_color = c("white", "#B20A2C"), opacity = 0.5,
                                                               #viridis::viridis(5)[1:4]),
                                                               #'#D7191C','#FDAE61','#FFFFBF','#A6D96A','#1A9641'),
                                                               #number_fmt = scales::number_format(scale=50), #scales::log_trans(), #scales::comma,
                                                               #bold_text = TRUE,
                                                               text_size = 10,
                                                               #background = '#555555',
                                                               min_value = 51, max_value=100,
                                                               #max_value = max(c(checkseq.sub$phred_raw,checkseq.sub$phred_trim)),
                                                               show_min_max = FALSE)),
                                # Barplot widgets
                                #-------------------
                                length = colDef(name = "trimLength", width=120, cell = data_bars(merged.df3,
                                                                                                    text_position = "inside-base",
                                                                                                    text_size = 10,
                                                                                                    #fill_color="#6f86ab",
                                                                                                    fill_color = c('#FFF2D9','#FFE1A6','#FFCB66','#FFB627'),
                                                                                                    fill_gradient = TRUE,
                                                                                                    background = 'transparent',
                                                                                                    number_fmt = scales::comma_format(accuracy = 0.1),
                                                                                                    round_edges = FALSE, align_bars="left"))
                             )) %>% add_title(paste("|assign_taxonomy| output for project: ", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))], sep=""))
                 
fname_html <- file.path(path, paste("output/assign_taxonomy_output_PASS___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))], ".html", sep=""))
suppressMessages(reactablefmtr::save_reactable_test(merged.df3_react, fname_html))
openFileInOS(fname_html)

#Clean up files
unlink(file.path(path,uc.out)) #Remove UC file
unlink(file.path(path,b6.out)) #Remove B6 file

#Export messages
message(cat(paste0("\033[97;", 40, "m","'assign_taxonomy' steps completed. Exporting files...", "\033[0m")))

message(cat(paste0("\033[97;", 40, "m","HTML file exported:", "\033[0m",
                   "\033[0;", 32, "m", " ", file.path(path, paste("output/assign_taxonomy_output_PASS___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))], ".html", sep="")),"\033[0m")))

if(export_csv==TRUE) {
  write.csv(merged.df3, file=paste("output/assign_taxonomy_output_PASS___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))],".csv", sep=""))
  message(cat(paste0("\033[97;", 40, "m","CSV file exported:", "\033[0m",
                     "\033[0;", 32, "m", " ", file.path(path, paste("output/assign_taxonomy_output_PASS___", unlist(strsplit(folder, '/'))[length(unlist(strsplit(folder, '/')))], ".csv", sep="")),"\033[0m", 
                     "\033[0;", 31, "m", "  <--- Input for Step 3: 'make_library'","\033[0m")))
}
}



#10^1.77
#}

#list(list(color_scales(merged.df3, color_ref = 'species_col', opacity = 0.5)), "white-space: nowrap", "font-size: 16px", "font-family: Calibri")
#lapply(merged.df3, function(x){
#       cell_style(merged.df3, font_size =10,  font_color = "white", rows = which(!is.na(merged.df3$species_col)), background_color = "#1673ba")
#       cell_style(merged.df3, font_size =10,  font_color = "white", rows = which(is.na(merged.df3$species_col)), background_color = "#1673ba")
#       })


#as.list(cell_style(merged.df3, font_size =10,  font_color = "white", rows = which(is.na(merged.df3$species_col)), background_color = "#1673ba"))

#style = cell_style(merged.df3, font_size =10,  font_color = "white", rows = which(is.na(merged.df3$species_col)), background_color = "#1673ba")),
#cell_style(merged.df3, font_size =10,  font_color = "white",
#            rows = which(is.na(merged.df3$species_col)), 
#          background_color = "grey"))),
#color_scales(merged.df3, color_ref = 'species_col', opacity = 0.5)),
#cell_style(merged.df3, font_size=20))), # "white-space: nowrap", "font-size: 16px", "font-family: Calibri")),
#style = as.list(color_scales(merged.df3, color_ref = 'species_col', opacity = 0.5), "white-space: nowrap", "font-size: 16px", "font-family: Calibri")),


#style = list(whiteSpace = "nowrap", fontFamily = "sans-serif")),
#style = list_with_names(color_scales(merged.df3, color_ref = 'species_colors', opacity = 0.5),
#                       whiteSpace = "nowrap", fontSize=10, fontFamily = "sans-serif")),
#style = list_with_names(color_scales(merged.df3, color_ref = 'species_colors', opacity = 0.5))),
