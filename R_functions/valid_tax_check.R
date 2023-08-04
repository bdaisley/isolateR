#valid_tax_check function

#Parameters:
# input = csv file (including full path if not in dir)
# col_species = column name with species names

#This function will determine if each species in a CSV file is validly published or not
#result file will be a csv with the results appended onto the input data

valid_tax_check<- function(input = NULL, col_species = "species"){
  if(require(LPSN)==FALSE) install.packages("LPSN", repos="http://R-Forge.R-project.org")
  if(require(stringr)==FALSE) install.packages("stringr")
  
  library("LPSN")
  library("stringr")
  
  if(is.null(input)) stop('Input file not supplied. Please check the execution script preamble to make sure it is entered correctly', call.=FALSE)
  
  #read file in
  inputdf <- read.csv(input, header = TRUE, stringsAsFactors = FALSE)
  
  #connect to the API
  message("To connect to the LPSN API an account a LPSN API account is required.")
  message("If you do not already possess one, please reigster here: https://api.lpsn.dsmz.de/")
  message("Please enter your LPSN API credentials")
  uname <- readline(prompt="Enter username: ")
  pword <- readline(prompt="Enter password: ")
  
  lpsn <- open_lpsn(uname, pword)
  stopifnot(
    inherits(lpsn, "lpsn_access"),
    !summary(lpsn)[c("expired", "refresh_expired")]
    ) #from fetch documentation, make sure session not expired
  
  
  #intialize dataframe for results
  checkdf <- as.data.frame(cbind(species = inputdf[,col_species], lpsnresult = NA, correctname = NA))
  
  #remove any tailing whitespace and pesky "Genus" or [Genus] naming conventions
  #also if user had _ seperating genus and species it will replace with space
  #if user had genus unclassified as species marker it will remove that and only check genus ID
  checkdf$species <- gsub("_", " ", checkdf$species)
  checkdf$species <- gsub("unclassified", "", checkdf$species)
  checkdf$species <- str_squish(checkdf$species)
  checkdf$species <- gsub("[[:punct:]]", "", checkdf$species)
  
  
  #for each species call to LPSN to retrieve info, if null not valid name
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
  
  outdf <- cbind(inputdf, LPSN_Result = checkdf$lpsnresult, Correct_Name = checkdf$correctname)
  
  outfile = paste0(gsub(".csv","",input),"_valid_species_check.csv")
  
  write.csv(outdf, outfile ,row.names = FALSE)
  
  
}
