#' @name make_fasta
#' @title Convert CSV file containing sequences to FASTA format
#' @description This function extracts sequences from a table in CSV format and converts them to FASTA format. Requires two columns, one with sequences and one with sequence names.
#' @export
#' @param csv_file Filename (or path and filename if not in working directory) of the table from which you would like to generate a FASTA file.
#' @param col_names Column name with the unique names/identifiers. (Default="ID")
#' @param col_seqs Column name with the sequences. (Default="Sequence")
#' @param output Desired filename for output FASTA file (Default = "output.fasta")
#' @return Returns sequences in FASTA format.
#' @importFrom BiocManager install
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings writeXStringSet
#' @importFrom utils read.csv
#' @importFrom stats setNames
#' @examples
#' #Set path to directory containing example .ab1 files
#' fpath1 <- system.file("extdata/abif_examples/rocket_salad", package = "isolateR")
#'
#' #Run isoQC function with default settings to generate CSV file
#' isoQC.S4 <- isoQC(input=fpath1)
#'
#' #Set path of CSV output file from isoQC step
#' csv.path <- file.path(fpath1, "isolateR_output/01_isoQC_trimmed_sequences_PASS.csv")
#'
#' #Run make_fasta function
#' make_fasta(csv_file= csv.path, col_names="filename", col_seqs="seqs_trim", output="output.fasta")

make_fasta <- function(csv_file=NULL,
                       col_names="ID",
                       col_seqs="Sequence",
                       output="output.fasta"){


  #load file----------------------
  df <- read.csv(csv_file, header = TRUE, stringsAsFactors = FALSE)

  #get data-------------------
  ss <- Biostrings::DNAStringSet(stats::setNames(df[,paste(col_seqs)], df[,paste(col_names)]))

  #output as fasta file-----------------
  fasta.output <- file.path(paste(unlist(strsplit(csv_file, '/'))[1:(length(unlist(strsplit(csv_file, '/')))-1)], collapse="/"), output)
  Biostrings::writeXStringSet(ss, fasta.output)
}
