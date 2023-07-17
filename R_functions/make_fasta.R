#make_fasta function

#Input--------------
#csv_file = filename (or path and filename if not in working dir)
#         of the table you would like to extract a fasta file from
#col_names = the column name with the unique names/identifiers
#col_seqs = the column name with the sequences



make_fasta <- function(csv_file=NULL,
                       col_names="ID",
                       col_seqs="Sequence"){
  
  #function requirements--------------------  
  if(require(Biostrings)==FALSE){
    if (!require("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("Biostrings", update = FALSE)
  }

  #load file----------------------
  df <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)
  
  
  #get data-------------------
  ss <- Biostrings::DNAStringSet(setNames(df[,paste(col_seqs)], df[,paste(col_names)]))
  
  #output as fasta file-----------------
  Biostrings::writeXStringSet(ss, "output.fasta" )
}