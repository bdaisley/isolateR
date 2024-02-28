#' @title Download taxonomic reference database
#' @name get_db
#' @description This function downloads taxonomic reference database and formats them for use.
#' @param db Database selection. One of "16S", "16S_arc", "18S", "ITS", or "cpn60"
#' @param force_update Forces new databases to be downloaded.
#' @returns Returns file path for database of interest
#' @importFrom utils download.file
#' @importFrom stringr str_subset
#' @importFrom R.utils gunzip
#' @examples
#' db.path <- get_db(db="16S", force_update=FALSE)

get_db <- function(db="16S", force_update=FALSE){

  #Check inputs----------------------------------------------------------------------------------------------------------------
  if(grepl("$16S_bac^|$16S_arc^|$18S_fun^|$cpn60^", db)) stop('Database entry incorrect. Please choose one of "16S_bac", "16S_arc", "18S_fun", "cpn60"')

  #Set paths-------------------------------------------------------------------------------------------------------------------
  db.path <- file.path(system.file("", package="isolateR"), "databases")
  suppressWarnings(dir.create(db.path))
  db.files  <- dir(db.path, full.names = FALSE)


  #FTP links for downloading reference DBs-------------------------------------------------------------------------------------
  dblink_16S <- "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
  dblink_16S_arc <- "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz"
  dblink_18S <- "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.18SrRNA.fna.gz"
  dblink_cpn60 <- "https://github.com/HillLabSask/cpn60-Classifier/releases/download/v11.1/cpn60-Classifier_v11.0_training.tar.gz"
  dblink_ITS <- "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz"

  #Set vector list for database matching---------------------------------------------------------------------------------------

  dblink.list <- setNames(c("dblink_16S", "dblink_16S_arc", "dblink_18S", "dblink_cpn60", "dblink_ITS"),
                          c("16S", "16S_arc", "18S", "cpn60", "ITS"))



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
  db.path <- file.path(db.path, paste(db, ".fna", sep=""))
  return(db.path)
}





