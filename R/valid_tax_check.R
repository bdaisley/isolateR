#' @name valid_tax_check
#' @title Validate species name via API client of LPSN
#' @description This function will determine if each species in a CSV file is validly published or not.
#' Result file will be a CSV with the results appended to the input data.
#' This function requires the user to have an LPSN API account setup. For more details and to register, see here: https://api.lpsn.dsmz.de/)
#' @export
#' @param input CSV file path. Expects full path if CSV file is not in the current working directory.
#' @param col_species Specify the column containing the binomial species names (e.g. "Akkermansia muciniphila")
#' @param export_csv Toggle (TRUE/FALSE). Set TRUE to automatically write .CSV file of results to current directory. (Default=TRUE)
#' @importFrom LPSN open_lpsn
#' @importFrom stringr str_squish
#' @return Returns a CSV saved in working directory

valid_tax_check<- function(input = NULL, col_species = "species", export_csv=TRUE){
  if(is.null(input)) stop('Input file not supplied. Please check the execution script preamble to make sure it is entered correctly', call.=FALSE)

  #read file in---------------------------------------------------------------------------------------------
  inputdf <- read.csv(input, header = TRUE, stringsAsFactors = FALSE)

  #connect to the API---------------------------------------------------------------------------------------
  message("To connect to the LPSN API an account a LPSN API account is required.")
  message("If you do not already possess one, please reigster here: https://api.lpsn.dsmz.de/")
  message("Please enter your LPSN API credentials")
  uname <- readline(prompt="Enter username: ")
  pword <- readline(prompt="Enter password: ")

  lpsn <- LPSN::open_lpsn(uname, pword)
  stopifnot(
    inherits(lpsn, "lpsn_access"),
    !summary(lpsn)[c("expired", "refresh_expired")]
    ) #from fetch documentation, make sure session not expired


  #intialize dataframe for results--------------------------------------------------------------------------
  checkdf <- as.data.frame(cbind(species = inputdf[,col_species], lpsnresult = NA, correctname = NA))

  #remove any tailing whitespace and pesky "Genus" or [Genus] naming conventions
  #also if user had _ seperating genus and species it will replace with space
  #if user had genus unclassified as species marker it will remove that and only check genus ID
  checkdf$species <- gsub("_", " ", checkdf$species)
  checkdf$species <- gsub("unclassified", "", checkdf$species)
  checkdf$species <- str_squish(checkdf$species)
  checkdf$species <- gsub("[[:punct:]]", "", checkdf$species)


  #for each species call to LPSN to retrieve info, if null not valid name----------------------------------
  for(i in 1:nrow(inputdf)){
    #use retrieve function to get information from full name
    ret <- retrieve(lpsn, search = "flexible", full_name = checkdf[i,"species"])

    #if there were no results, set as non valid and move on
    if(length(ret) == 0){
      checkdf[i, "lpsnresult"] <- "Not valid"
      next
    }

    #if there is a result, check if it is a synonym
    retdf <- as.data.frame(ret)
    if(retdf[1,"lpsn_taxonomic_status"] == "synonym"){
      retsyn <- as.data.frame(retrieve(lpsn, search = "flexible", id = retdf[1,"lpsn_correct_name_id"]))
      checkdf[i, "lpsnresult"] <- "Synonym"
      checkdf[i, "correctname"] <- retsyn[1,"full_name"]
      next
    }

    if(retdf[1,"validly_published"] == "ICNP"){
      checkdf[i, "lpsnresult"] <- "Validly published"
      checkdf[i, "correctname"] <- retdf[1,"full_name"]
    } else {
      #just to catch anything weird going on, nothing should get to this point
      message(paste("Something weird with this species ",checkdf[i,"species"]))
      checkdf[i, "lpsnresult"] <- retdf[1,"validly_published"]
      checkdf[i, "correctname"] <- retdf[1,"full_name"]
    }

  }

  #output------------------------------------------------------------------------------------------------
  outdf <- cbind(inputdf, LPSN_Result = checkdf$lpsnresult, Correct_Name = checkdf$correctname)

  outfile = paste0(gsub(".csv","",input),"_valid_species_check.csv")

  if(export_csv==TRUE){
    write.csv(outdf, outfile ,row.names = FALSE)
  }

  return(outdf)

}
