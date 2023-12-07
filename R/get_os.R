#' @name get_os
#' @title Determine user operating system.
#' @description Determines the type of operating system being used.
#' @export
#' @return Returns sysname as one of windows/osx-mac/linux
#' @examples
#' #Example 1 on a Windows-based operating system
#' os.index <- get_os()
#' print(os.index)
#'
#' #Example 2 on a Mac operating system
#' os.index <- get_os()
#' print(os.index)
#'
#' #Example 3 on a Linux operating system
#' os.index <- get_os()
#' print(os.index)


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
