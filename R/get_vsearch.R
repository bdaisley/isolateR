#' @name get_vsearch
#' @title Download VSEARCH software reference database
#' @description This function downloads the VSEARCH software used querying sequences against taxonomic databases of interest.
#' @param os Operating system, one of: "windows", "osx-mac", or "linux". If blank (os=NULL) then will try to automatically determine operating system.
#' @export
#' @returns Returns path for VSEARCH executable
#' @importFrom utils download.file
#' @importFrom utils unzip
#' @importFrom stringr str_subset
#' @importFrom R.utils gunzip
#' @examples
#' #Example for automatically detecting operating system and downloading VSEARCH software
#' vsearch.path <- get_vsearch()

get_vsearch <- function(os=NULL){
  message(cat(paste0("\n", "\033[97;", 40, "m","Detecting operating system...", "\033[0m")))
  
  vsearch.path.dl <- file.path(system.file(package="isolateR"), "vsearch")
  suppressWarnings(dir.create(vsearch.path.dl))
  vsearch_files <- stringr::str_subset(dir(vsearch.path.dl, full.names = FALSE), 'vsearch')

  #---Get operating system
  if(is.null(os)){os <- paste(isolateR::get_os())}
  
  if(identical(vsearch_files, character(0))){
    if(os=="windows"){
      if(identical(vsearch_files, character(0))){
        message(cat(paste0("\n", "\033[0;", 32, "m","Operating system is ---> Windows-based <---", "\033[0m")))
        utils::download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-win-x86_64.zip", file.path(vsearch.path.dl, 'vsearch-2.23.0-win-x86_64.zip'), mode='wb')
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
    
    if(os=="osx-mac"){
      if(identical(vsearch_files, character(0))){
        message(cat(paste0("\033[0;", 32, "m","Operating system is ---> MacOS-based <---", "\033[0m")))
        utils::download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-macos-universal.tar.gz", file.path(vsearch.path.dl, 'vsearch-2.23.0-macos-universal.tar.gz'), mode='wb')
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
    
    if(os=="linux"){
      if(identical(vsearch_files, character(0))){
        message(cat(paste0("\033[0;", 32, "m","Operating system is ---> Linux-based <---", "\033[0m")))
        utils::download.file("https://github.com/torognes/vsearch/releases/download/v2.23.0/vsearch-2.23.0-linux-x86_64.tar.gz", file.path(vsearch.path.dl, 'vsearch-2.23.0-linux-x86_64.tar.gz'), mode='wb')
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
  } else{
    message(cat(paste0("\033[0;", 32, "m","VSEARCH already downloaded.", "\033[0m")))
    if(os=="windows"){vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0.exe")}
    if(os=="osx-mac"){vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0_macos")}
    if(os=="linux"){vsearch.path <- file.path(vsearch.path.dl,"vsearch-2.23.0")}
  }
  return(vsearch.path)
}
