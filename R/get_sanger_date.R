#' @name get_sanger_date
#' @title get_sanger_date function
#' @description Helper function to automatically retrieve run date from Sanger sequencing .ab1 files.
#' @export
#' @param file The .ab1 file in from which to retrieve the date information. (Must be in S4 abif format)
#' @return Returns date in "YYYY_MM_DD" format
#' @importFrom stringr str_pad
#' @importFrom sangerseqR read.abif
#' @examples
#' #Path to the first listed .ab1 file in example directory
#' fpath <- file.path(system.file("extdata/abif_examples/rocket_salad", package = "isolateR"),
#'                    list.files(system.file("extdata/abif_examples/rocket_salad", package = "isolateR"))[1])
#' #Read in the ab1 file to S4 format
#' ab1.S4 <- sangerseqR::read.abif(fpath)
#'
#' #Retrieve date
#' get_sanger_date(ab1.S4)


get_sanger_date <- function(file=NULL){
sanger.date <- c(paste(file@data['RUND.1'][[1]][[1]],
                       stringr::str_pad(file@data['RUND.1'][[1]][[2]], 2, pad = "0"),
                       stringr::str_pad(file@data['RUND.1'][[1]][[3]], 2, pad = "0"), sep="_"))
return(sanger.date)
}
