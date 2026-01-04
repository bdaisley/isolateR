#' @title Download taxonomic reference database
#' @name get_db
#' @description This function downloads taxonomic reference database and formats them for use.
#' @param db Database selection. One of "16S", "16S_arc", "18S", "ITS", or "cpn60"
#' @param force_update Forces new databases to be downloaded.
#' @param add_taxonomy Add full taxonomy to header to enable classification of sequences offline. This is primarily for cases of operation with no internet or when computing within an HPC cluster without administrator access. Results in a semicolon (;) delimited FASTA header as follows: >Accession_no;d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species (e.g., >NR_042817.1;d__Bacteria;p__Verrucomicrobiota;c__Verrucomicrobiia;o__Verrucomicrobiales;f__Akkermansiaceae;g__Akkermansia;s__Akkermansia_muciniphila)
#' @export
#' @returns Returns file path for database of interest
#' @importFrom utils download.file
#' @importFrom stringr str_subset
#' @importFrom R.utils gunzip
#' @examples
#' db.path <- get_db(db="16S", force_update=FALSE)

get_db <- function(db="16S_bac", force_update=FALSE, add_taxonomy=FALSE){

  #Check inputs----------------------------------------------------------------------------------------------------------------
  if(!grepl("^16S$|^16S_bac$|^16S_arc$|^18S_fun$|^cpn60$|^ITS$", db)) stop('Database entry incorrect. Please choose one of: "16S_bac", "16S_arc", "18S_fun", "cpn60", "ITS"')
  if(grepl("^16S$", db)) message('Note: Use of db="16S" defaults to the Targeted Loci Bacteria 16S rRNA database (db="16S_bac"). Please use db="16_arc" for Archaea.')
  #Set paths-------------------------------------------------------------------------------------------------------------------
  db.path <- file.path(system.file("", package="isolateR"), "databases")
  suppressWarnings(dir.create(db.path))
  db.files  <- dir(db.path, full.names = FALSE)
  
  
  #FTP links for downloading reference DBs-------------------------------------------------------------------------------------
  dblink_16S_bac <- "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
  dblink_16S_arc <- "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz"
  dblink_18S_fun <- "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.18SrRNA.fna.gz"
  dblink_cpn60 <- "https://github.com/HillLabSask/cpn60-Classifier/releases/download/v11.1/cpn60-Classifier_v11.0_training.tar.gz"
  dblink_ITS <- "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz"
  
  #Set vector list for database matching---------------------------------------------------------------------------------------
  
  dblink.list <- setNames(c("dblink_16S_bac", "dblink_16S_bac", "dblink_16S_arc", "dblink_18S_fun", "dblink_18S_fun", "dblink_cpn60", "dblink_ITS"),
                          c("16S_bac", "16S", "16S_arc", "18S", "18S_fun", "cpn60", "ITS"))
  
  
  
  #cpn60 site php specifics----------------------------------------------------------------------------------------------------
  head.list <- NULL
  if(grepl("cpn60", db)){
    head.list <- setNames(c("www.cpndb.ca",
                            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/119.0.0.0 Safari/537.36",
                            "text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,image/apng,*/*;q=0.8,application/signed-exchange;v=b3;q=0.7",
                            "en-US,en;q=0.9",
                            "https://www.cpndb.ca/downloads.php",
                            "PHPSESSID=1fuqm20nk0fmo6auiujbrhdf7i",
                            "keep-alive"),
                          c("Host",
                            "User-Agent",
                            "Accept",
                            "Accept-Language",
                            "Referer",
                            "Cookie",
                            "Connection"))
  }
  
  
  #Download db---------------------------------------------------------------------------------------
  
  if(force_update==FALSE){
    if(!(isTRUE(paste(db, ".fna", sep="") == (stringr::str_subset(dir(db.path, full.names = FALSE), paste(db, ".fna", sep="")))) )){
      if(db == "cpn60"){
        utils::download.file(eval(parse(text = paste(dblink.list)[names(dblink.list) %in% db])),
                             file.path(db.path, paste(db, ".tar.gz", sep="")),
                             method= "libcurl",
                             mode='wb',
                             headers=NULL)
        R.utils::gunzip(file.path(db.path, paste(db, ".tar.gz", sep="")), remove = FALSE, overwrite=TRUE)
        utils::untar(file.path(db.path, paste(db, ".tar", sep="")), files= "cpn60-Classifier_v11.0_training/refseqs_v11.fasta" , exdir=file.path(db.path))
        file.rename(file.path(db.path, paste("cpn60-Classifier_v11.0_training/refseqs_v11.fasta", sep="")), file.path(db.path, paste(db, ".fna", sep="")))
        unlink(file.path(db.path, paste("cpn60-Classifier_v11.0_training", sep="")) ,recursive=TRUE)
        unlink(file.path(db.path, paste(db, ".tar.gz", sep="")), recursive=TRUE)
        unlink(file.path(db.path, paste(db, ".tar", sep="")), recursive=TRUE)
        db.fasta.in <- Biostrings::readBStringSet(file.path(db.path, paste(db, ".fna", sep="")))
        names(db.fasta.in) <- stringr::str_split_fixed(names(db.fasta.in), " ", 4)[,2]
        Biostrings::writeXStringSet(db.fasta.in, filepath=file.path(db.path, paste(db, ".fna", sep="")) , format="fasta")
      } else {
        utils::download.file(eval(parse(text = paste(dblink.list)[names(dblink.list) %in% db])),
                             file.path(db.path, paste(db, ".fna.gz", sep="")),
                             method= "libcurl",
                             mode='wb',
                             headers=NULL)
        R.utils::gunzip(file.path(db.path, paste(db, ".fna.gz", sep="")), remove = FALSE, overwrite=TRUE)
        db.fasta.path <- file.path(db.path, paste(db, ".fna", sep=""))
        unlink(file.path(db.path, paste(db, ".fna.gz", sep="")),recursive=TRUE)
      }
    }
  }
  
  if(force_update==TRUE){
    if(db == "cpn60"){
      utils::download.file(eval(parse(text = paste(dblink.list)[names(dblink.list) %in% db])),
                           file.path(db.path, paste(db, ".tar.gz", sep="")),
                           method= "libcurl",
                           mode='wb',
                           headers=NULL)
      R.utils::gunzip(file.path(db.path, paste(db, ".tar.gz", sep="")), remove = FALSE, overwrite=TRUE)
      utils::untar(file.path(db.path, paste(db, ".tar", sep="")), files= "cpn60-Classifier_v11.0_training/refseqs_v11.fasta" , exdir=file.path(db.path))
      file.rename(file.path(db.path, paste("cpn60-Classifier_v11.0_training/refseqs_v11.fasta", sep="")), file.path(db.path, paste(db, ".fna", sep="")))
      unlink(file.path(db.path, paste("cpn60-Classifier_v11.0_training", sep="")) ,recursive=TRUE)
      unlink(file.path(db.path, paste(db, ".tar.gz", sep="")), recursive=TRUE)
      unlink(file.path(db.path, paste(db, ".tar", sep="")), recursive=TRUE)
      db.fasta.in <- Biostrings::readBStringSet(file.path(db.path, paste(db, ".fna", sep="")))
      names(db.fasta.in) <- stringr::str_split_fixed(names(db.fasta.in), " ", 4)[,2]
      Biostrings::writeXStringSet(db.fasta.in, filepath=file.path(db.path, paste(db, ".fna", sep="")) , format="fasta")
    } else {
      utils::download.file(eval(parse(text = paste(dblink.list)[names(dblink.list) %in% db])),
                           file.path(db.path, paste(db, ".fna.gz", sep="")),
                           method= "libcurl",
                           mode='wb',
                           headers=NULL)
      R.utils::gunzip(file.path(db.path, paste(db, ".fna.gz", sep="")), remove = FALSE, overwrite=TRUE)
      db.fasta.path <- file.path(db.path, paste(db, ".fna", sep=""))
      unlink(file.path(db.path, paste(db, ".fna.gz", sep="")),recursive=TRUE)
    }
  }
  
  db.fasta.path <- file.path(db.path, paste(db, ".fna", sep=""))
  
  #Add taxonomy---------------------------------------------------------------------------------------
  if(force_update==FALSE){
    if(add_taxonomy==TRUE){
      if(!(isTRUE(paste(db, "___taxonomy_in_headers.fna", sep="") == (stringr::str_subset(dir(db.path, full.names = FALSE), paste(db, "___taxonomy_in_headers.fna", sep="")))) )){
        db.fasta <- Biostrings::readBStringSet(db.fasta.path)
        acc.list <- stringr::str_split_fixed(names(db.fasta), " ", 2)[,1]
        #:::::::::::::::::::::::::::::::::::::
        # Fetch NCBI metadata for sequences
        #:::::::::::::::::::::::::::::::::::::
        id.list <- 1:length(acc.list)
        id.chunk <- 250
        fetch.ids <- split(id.list, ceiling(seq_along(id.list)/id.chunk))
        
        message(cat(paste0("\033[97;", 40, "m","Downloading full taxonomy for ", length(acc.list), " sequences. This may take a while...", "\033[0m", "\n")))
        fetch.list <- list()
        for (i in 1:length(fetch.ids)) {
          ii <- fetch.ids[[i]]
          while(TRUE){
            links <- rentrez::entrez_link(dbfrom="nucleotide", id=acc.list[ii], db="taxonomy")
            taxids <- links$links$nuccore_taxonomy
            r_fetch.try <- try(rentrez::entrez_fetch(db="taxonomy", id=taxids, rettype="xml")) # OR use rettype="xml" OR useing rettype="gbc", silent=TRUE)
            if(is(r_fetch.try, 'try-error'))
              cat("Failed, trying again in 10 seconds...\n")
            Sys.sleep(3)
            if(!is(r_fetch.try, 'try-error')) break
          }
          while(TRUE){
            r_fetch.try2 <- try(rentrez::entrez_fetch(db="nucleotide", id=acc.list[ii], rettype="gbc")) # OR use rettype="xml" OR useing rettype="gbc", silent=TRUE)
            if(is(r_fetch.try2, 'try-error'))
              cat("Failed, trying again in 10 seconds...\n")
            Sys.sleep(3)
            if(!is(r_fetch.try2, 'try-error')) break
          }
          
          #---------- Get type strain culture collection numbers
          r_fetch2 <- stringr::str_replace_all(r_fetch.try2, "INSDSeq_update-date", "INSDSeq_update_date")
          r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_create-date", "INSDSeq_create_date")
          r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_primary.accession", "INSDSeq_primary_accession")
          r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_accession.version", "INSDSeq_accession_version")
          r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_other.seqids", "INSDSeq_other_seqids")
          r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_feature-table", "INSDSeq_feature_table")
          r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_feature-table", "INSDSeq_feature_table")
          
          r_fetch2.df <- xmlconvert::xml_to_df(text = r_fetch2, records.tags = c("//INSDSeq"), check.datatypes=FALSE, check.names=TRUE) %>%
            dplyr::rowwise() %>%
            mutate(taxid = unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'INSDQualifier_value~taxon:', fixed=TRUE))[2], '|', fixed=TRUE))[1]) %>%
            mutate(strain = unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'strain||INSDQualifier_value~', fixed=TRUE))[2], '|', fixed=TRUE))[1]) %>%
            mutate(culture_collection = gsub(":", " ", unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'culture_collection||INSDQualifier_value~', fixed=TRUE))[2], '|', fixed=TRUE))[1])) %>%
            mutate(culture_collection = ifelse(is.na(culture_collection), strain, culture_collection)) %>%
            ungroup() %>%
            mutate(culture_collection = ifelse(is.na(strain) & is.na(culture_collection), 
                                               stringr::str_split_fixed(stringr::str_extract(INSDSeq_definition, " [A-Z]+.*"), " ITS | 16S ", 2)[,1],
                                               culture_collection)) %>%
            mutate(culture_collection = gsub("^ ", "", culture_collection)) %>%
            mutate(closest_match = paste(stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,1], " ", stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,2], " (", culture_collection, ")", sep=""))
          
          #---------- #Get taxonomy lineage info
          r_fetch <- r_fetch.try
          r_fetch.xml <- xml2::read_xml(r_fetch)
          top_taxa <- xml2::xml_find_all(r_fetch.xml, "//TaxaSet/Taxon")
          
          all_lineages <- lapply(1:length(top_taxa), function(taxon_node){
            lineage_nodes <- xml2::xml_find_all(top_taxa[taxon_node], "./LineageEx/Taxon")
            #------------------------------------------------------------------
            main_taxid <- xml2::xml_text(xml2::xml_find_all(top_taxa[taxon_node], "./TaxId"))
            main_species <- xml2::xml_text(xml2::xml_find_all(top_taxa[taxon_node], "./ScientificName"))
            main_rank <- xml2::xml_text(xml2::xml_find_all(top_taxa[taxon_node], "./Rank"))
            #------------------------------------------------------------------
            taxname = xml2::xml_text(xml2::xml_find_all(lineage_nodes, "./ScientificName"))
            taxrank = xml2::xml_text(xml2::xml_find_all(lineage_nodes, "./Rank"))
            #------------------------------------------------------------------
            tibble(taxid = main_taxid,
                   #INSDSeq_accession_version = r_fetch2.df$INSDSeq_accession_version[taxon_node],
                   #accession = r_fetch2.df$INSDSeq_accession_version[taxon_node],
                   #closest_match = r_fetch2.df$closest_match[taxon_node],
                   domain =  ifelse(isEmpty(taxname[taxrank == "domain"]), "NA", taxname[taxrank == "domain"]),
                   phylum =  ifelse(isEmpty(taxname[taxrank == "phylum"]), "NA", taxname[taxrank == "phylum"]),
                   class =  ifelse(isEmpty(taxname[taxrank == "class"]), "NA", taxname[taxrank == "class"]),
                   order = ifelse(isEmpty(taxname[taxrank == "order"]), "NA", taxname[taxrank == "order"]),
                   family = ifelse(isEmpty(taxname[taxrank == "family"]), "NA", taxname[taxrank == "family"]),
                   genus = ifelse(isEmpty(taxname[taxrank == "genus"]), "NA", taxname[taxrank == "genus"]),
                   species = ifelse(isEmpty(main_species), "NA", main_species))
            
          })
          
          all_lineages.df <- dplyr::bind_rows(all_lineages, .id = "column_label") %>% filter(taxid!=1) #%>% dplyr::relocate("accession", 1)
          if(i == 1){
            fetch.list <- all_lineages.df
            fetch.list2 <- r_fetch2.df}
          if(i > 1){
            fetch.list <- rbind(fetch.list, all_lineages.df)
            fetch.list2 <- rbind(fetch.list2, r_fetch2.df)
          }
          svMisc::progress(max(fetch.ids[[i]]), length(acc.list))
        }
        
        #Add columns of interest to lookup table
        suppressWarnings(fetch.list.df <-  merge(fetch.list2, fetch.list, by="taxid", sort=FALSE, all=TRUE) %>%
                           #mutate(closest_match = paste(stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,1], " ", stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,2], " (", culture_collection, ")", sep="")) %>%
                           mutate(species = paste(stringr::str_split_fixed(closest_match, " ", 3)[,1],stringr::str_split_fixed(closest_match, " ", 3)[,2], sep=" ")) %>%
                           mutate(taxonomy = paste(domain, phylum, class, order, family, genus, species, sep=";")) %>%
                           mutate(rank_domain = domain) %>%
                           mutate(rank_phylum = phylum) %>%
                           mutate(rank_class = class) %>%
                           mutate(rank_order = order) %>%
                           mutate(rank_family = family) %>%
                           mutate(rank_genus = genus) %>%
                           mutate(rank_species = gsub("'", "", species)) %>% #Replace instances where single quotations are in species name
                           #mutate(genus_tmp = stringr::str_split_fixed(rank_species, " ", 2)[,1]) %>% #Fixing instances where genus is in wrong spot
                           #mutate(genus_tmp = gsub("\\[|\\]", "", genus_tmp)) %>% #Fixing instances where genus is in wrong spot
                           mutate(species = gsub(" ", "_", species)) %>% #Replace spaces in species name
                           mutate(rank_genus = stringr::str_split_fixed(species, "_", 2)[,1]) %>%
                           mutate(INSDSeq_taxonomy = paste(rank_domain, rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, sep=";")) %>%
                           mutate_at(vars(rank_class, rank_order, rank_family), funs(ifelse(. == "", "NA", .))) #Replace unknown ranks with "NA"
        )
        
        lookup.list <- setNames(fetch.list.df$INSDSeq_taxonomy, fetch.list.df$INSDSeq_accession_version)
        lookup.list2 <- setNames(fetch.list.df$closest_match, fetch.list.df$INSDSeq_accession_version)
        lookup.list3 <- setNames(fetch.list.df$species, fetch.list.df$INSDSeq_accession_version)

        fetch.list.df.mer.x <- as.data.frame(cbind("NCBI_acc" = acc.list)) %>% 
          mutate(taxonomy = dplyr::recode(NCBI_acc, !!!lookup.list)) %>%
          mutate(taxonomy = gsub(" ", "", .$taxonomy)) %>%
          mutate(rank_domain = stringr::str_split_fixed(.$taxonomy, ";", 7)[,1]) %>%
          mutate(rank_phylum = stringr::str_split_fixed(.$taxonomy, ";", 7)[,2]) %>%
          mutate(rank_class = stringr::str_split_fixed(.$taxonomy, ";", 7)[,3]) %>%
          mutate(rank_order = stringr::str_split_fixed(.$taxonomy, ";", 7)[,4]) %>%
          mutate(rank_family = stringr::str_split_fixed(.$taxonomy, ";", 7)[,5]) %>%
          mutate(rank_genus = stringr::str_split_fixed(.$taxonomy, ";", 7)[,6]) %>%
          mutate(species = dplyr::recode(NCBI_acc, !!!lookup.list3))
        
        #Add taxonomy to fasta header for database sequences
        names(db.fasta) <- paste(fetch.list.df.mer.x$NCBI_acc,
                                 ";d__", fetch.list.df.mer.x$rank_domain, 
                                 ";p__", fetch.list.df.mer.x$rank_phylum, 
                                 ";c__", fetch.list.df.mer.x$rank_class, 
                                 ";o__", fetch.list.df.mer.x$rank_order, 
                                 ";f__", fetch.list.df.mer.x$rank_family, 
                                 ";g__", fetch.list.df.mer.x$rank_genus, 
                                 ";s__", fetch.list.df.mer.x$species, sep="")
        
        #Write database to FASTA file
        db.fasta.path <- file.path(db.path, paste(db, "___taxonomy_in_headers.fna", sep=""))
        Biostrings::writeXStringSet(db.fasta, file=db.fasta.path)
      }
      #Return path if file already present
      db.fasta.path <- file.path(db.path, paste(db, "___taxonomy_in_headers.fna", sep=""))
    }
  }
  
  if(force_update==TRUE){
    if(add_taxonomy==TRUE){
      db.fasta <- Biostrings::readBStringSet(db.fasta.path)
      acc.list <- stringr::str_split_fixed(names(db.fasta), " ", 2)[,1]
      #:::::::::::::::::::::::::::::::::::::
      # Fetch NCBI metadata for sequences
      #:::::::::::::::::::::::::::::::::::::
      id.list <- 1:length(acc.list)
      id.chunk <- 250
      fetch.ids <- split(id.list, ceiling(seq_along(id.list)/id.chunk))
      
      message(cat(paste0("\033[97;", 40, "m","Downloading full taxonomy for ", length(acc.list), " sequences. This may take a while...", "\033[0m", "\n")))
      fetch.list <- list()
      for (i in 1:length(fetch.ids)) {
        ii <- fetch.ids[[i]]
        while(TRUE){
          links <- rentrez::entrez_link(dbfrom="nucleotide", id=acc.list[ii], db="taxonomy")
          taxids <- links$links$nuccore_taxonomy
          r_fetch.try <- try(rentrez::entrez_fetch(db="taxonomy", id=taxids, rettype="xml")) # OR use rettype="xml" OR useing rettype="gbc", silent=TRUE)
          if(is(r_fetch.try, 'try-error'))
            cat("Failed, trying again in 10 seconds...\n")
          Sys.sleep(3)
          if(!is(r_fetch.try, 'try-error')) break
        }
        while(TRUE){
          r_fetch.try2 <- try(rentrez::entrez_fetch(db="nucleotide", id=acc.list[ii], rettype="gbc")) # OR use rettype="xml" OR useing rettype="gbc", silent=TRUE)
          if(is(r_fetch.try2, 'try-error'))
            cat("Failed, trying again in 10 seconds...\n")
          Sys.sleep(3)
          if(!is(r_fetch.try2, 'try-error')) break
        }
        
        #---------- Get type strain culture collection numbers
        r_fetch2 <- stringr::str_replace_all(r_fetch.try2, "INSDSeq_update-date", "INSDSeq_update_date")
        r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_create-date", "INSDSeq_create_date")
        r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_primary.accession", "INSDSeq_primary_accession")
        r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_accession.version", "INSDSeq_accession_version")
        r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_other.seqids", "INSDSeq_other_seqids")
        r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_feature-table", "INSDSeq_feature_table")
        r_fetch2 <- stringr::str_replace_all(r_fetch2, "INSDSeq_feature-table", "INSDSeq_feature_table")
        
        r_fetch2.df <- xmlconvert::xml_to_df(text = r_fetch2, records.tags = c("//INSDSeq"), check.datatypes=FALSE, check.names=TRUE) %>%
          dplyr::rowwise() %>%
          mutate(taxid = unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'INSDQualifier_value~taxon:', fixed=TRUE))[2], '|', fixed=TRUE))[1]) %>%
          mutate(strain = unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'strain||INSDQualifier_value~', fixed=TRUE))[2], '|', fixed=TRUE))[1]) %>%
          mutate(culture_collection = gsub(":", " ", unlist(strsplit(unlist(strsplit(INSDSeq_feature_table, 'culture_collection||INSDQualifier_value~', fixed=TRUE))[2], '|', fixed=TRUE))[1])) %>%
          mutate(culture_collection = ifelse(is.na(culture_collection), strain, culture_collection)) %>%
          ungroup() %>%
          mutate(culture_collection = ifelse(is.na(strain) & is.na(culture_collection), 
                                             stringr::str_split_fixed(stringr::str_extract(INSDSeq_definition, " [A-Z]+.*"), " ITS | 16S ", 2)[,1],
                                             culture_collection)) %>%
          mutate(culture_collection = gsub("^ ", "", culture_collection)) %>%
          mutate(closest_match = paste(stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,1], " ", stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,2], " (", culture_collection, ")", sep=""))
        
        #---------- #Get taxonomy lineage info
        r_fetch <- r_fetch.try
        r_fetch.xml <- xml2::read_xml(r_fetch)
        top_taxa <- xml2::xml_find_all(r_fetch.xml, "//TaxaSet/Taxon")
        
        all_lineages <- lapply(1:length(top_taxa), function(taxon_node){
          lineage_nodes <- xml2::xml_find_all(top_taxa[taxon_node], "./LineageEx/Taxon")
          #------------------------------------------------------------------
          main_taxid <- xml2::xml_text(xml2::xml_find_all(top_taxa[taxon_node], "./TaxId"))
          main_species <- xml2::xml_text(xml2::xml_find_all(top_taxa[taxon_node], "./ScientificName"))
          main_rank <- xml2::xml_text(xml2::xml_find_all(top_taxa[taxon_node], "./Rank"))
          #------------------------------------------------------------------
          taxname = xml2::xml_text(xml2::xml_find_all(lineage_nodes, "./ScientificName"))
          taxrank = xml2::xml_text(xml2::xml_find_all(lineage_nodes, "./Rank"))
          #------------------------------------------------------------------
          tibble(taxid = main_taxid,
                 #INSDSeq_accession_version = r_fetch2.df$INSDSeq_accession_version[taxon_node],
                 #accession = r_fetch2.df$INSDSeq_accession_version[taxon_node],
                 #closest_match = r_fetch2.df$closest_match[taxon_node],
                 domain =  ifelse(isEmpty(taxname[taxrank == "domain"]), "NA", taxname[taxrank == "domain"]),
                 phylum =  ifelse(isEmpty(taxname[taxrank == "phylum"]), "NA", taxname[taxrank == "phylum"]),
                 class =  ifelse(isEmpty(taxname[taxrank == "class"]), "NA", taxname[taxrank == "class"]),
                 order = ifelse(isEmpty(taxname[taxrank == "order"]), "NA", taxname[taxrank == "order"]),
                 family = ifelse(isEmpty(taxname[taxrank == "family"]), "NA", taxname[taxrank == "family"]),
                 genus = ifelse(isEmpty(taxname[taxrank == "genus"]), "NA", taxname[taxrank == "genus"]),
                 species = ifelse(isEmpty(main_species), "NA", main_species))
          
        })
        
        all_lineages.df <- dplyr::bind_rows(all_lineages, .id = "column_label") %>% filter(taxid!=1) #%>% dplyr::relocate("accession", 1)
        if(i == 1){
          fetch.list <- all_lineages.df
          fetch.list2 <- r_fetch2.df}
        if(i > 1){
          fetch.list <- rbind(fetch.list, all_lineages.df)
          fetch.list2 <- rbind(fetch.list2, r_fetch2.df)
        }
        svMisc::progress(max(fetch.ids[[i]]), length(acc.list))
      }
      
      #Add columns of interest to lookup table
      suppressWarnings(fetch.list.df <-  merge(fetch.list2, fetch.list, by="taxid", sort=FALSE, all=TRUE) %>%
                         #mutate(closest_match = paste(stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,1], " ", stringr::str_split_fixed(INSDSeq_organism, " ", 3)[,2], " (", culture_collection, ")", sep="")) %>%
                         mutate(species = paste(stringr::str_split_fixed(closest_match, " ", 3)[,1],stringr::str_split_fixed(closest_match, " ", 3)[,2], sep=" ")) %>%
                         mutate(taxonomy = paste(domain, phylum, class, order, family, genus, species, sep=";")) %>%
                         mutate(rank_domain = domain) %>%
                         mutate(rank_phylum = phylum) %>%
                         mutate(rank_class = class) %>%
                         mutate(rank_order = order) %>%
                         mutate(rank_family = family) %>%
                         mutate(rank_genus = genus) %>%
                         mutate(rank_species = gsub("'", "", species)) %>% #Replace instances where single quotations are in species name
                         #mutate(genus_tmp = stringr::str_split_fixed(rank_species, " ", 2)[,1]) %>% #Fixing instances where genus is in wrong spot
                         #mutate(genus_tmp = gsub("\\[|\\]", "", genus_tmp)) %>% #Fixing instances where genus is in wrong spot
                         mutate(species = gsub(" ", "_", species)) %>% #Replace spaces in species name
                         mutate(rank_genus = stringr::str_split_fixed(species, "_", 2)[,1]) %>%
                         mutate(INSDSeq_taxonomy = paste(rank_domain, rank_phylum, rank_class, rank_order, rank_family, rank_genus, rank_species, sep=";")) %>%
                         mutate_at(vars(rank_class, rank_order, rank_family), funs(ifelse(. == "", "NA", .))) #Replace unknown ranks with "NA"
      )
      
      lookup.list <- setNames(fetch.list.df$INSDSeq_taxonomy, fetch.list.df$INSDSeq_accession_version)
      lookup.list2 <- setNames(fetch.list.df$closest_match, fetch.list.df$INSDSeq_accession_version)
      lookup.list3 <- setNames(fetch.list.df$species, fetch.list.df$INSDSeq_accession_version)
      
      fetch.list.df.mer.x <- as.data.frame(cbind("NCBI_acc" = acc.list)) %>% 
        mutate(taxonomy = dplyr::recode(NCBI_acc, !!!lookup.list)) %>%
        mutate(taxonomy = gsub(" ", "", .$taxonomy)) %>%
        mutate(rank_domain = stringr::str_split_fixed(.$taxonomy, ";", 7)[,1]) %>%
        mutate(rank_phylum = stringr::str_split_fixed(.$taxonomy, ";", 7)[,2]) %>%
        mutate(rank_class = stringr::str_split_fixed(.$taxonomy, ";", 7)[,3]) %>%
        mutate(rank_order = stringr::str_split_fixed(.$taxonomy, ";", 7)[,4]) %>%
        mutate(rank_family = stringr::str_split_fixed(.$taxonomy, ";", 7)[,5]) %>%
        mutate(rank_genus = stringr::str_split_fixed(.$taxonomy, ";", 7)[,6]) %>%
        mutate(species = dplyr::recode(NCBI_acc, !!!lookup.list3))
      
      #Add taxonomy to fasta header for database sequences
      names(db.fasta) <- paste(fetch.list.df.mer.x$NCBI_acc,
                               ";d__", fetch.list.df.mer.x$rank_domain, 
                               ";p__", fetch.list.df.mer.x$rank_phylum, 
                               ";c__", fetch.list.df.mer.x$rank_class, 
                               ";o__", fetch.list.df.mer.x$rank_order, 
                               ";f__", fetch.list.df.mer.x$rank_family, 
                               ";g__", fetch.list.df.mer.x$rank_genus, 
                               ";s__", fetch.list.df.mer.x$species, sep="")
      
      #Write database to FASTA file
      db.fasta.path <- file.path(db.path, paste(db, "___taxonomy_in_headers.fna", sep=""))
      Biostrings::writeXStringSet(db.fasta, file=db.fasta.path)
    }
    #Return path if file already present
    db.fasta.path <- file.path(db.path, paste(db, "___taxonomy_in_headers.fna", sep=""))
  }
  
  #message("Database path:")
  return(db.fasta.path)
}
