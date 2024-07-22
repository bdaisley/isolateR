###########################################################################################################

#::::::::::::::::::::::::::::::::::::::::
# Set 'class-isoQC' Class Object
#::::::::::::::::::::::::::::::::::::::::
#' @name class-isoQC
#' @rdname class-isoQC
#' @title isoQC Class Object
#' @description S4 wrapper for \code{\link{isoQC}} function. Access data via S4 slot functions.
#' @export
#' @slot date Character string containing run date from each of the input Sanger sequence .ab1 files ("YYYY_MM_DD" format).
#' @slot filename Character string containing input filenames.
#' @slot trim.start.pos Numeric string containing trimming position start point.
#' @slot trim.end.pos Numeric string containing trimming position end point.
#' @slot phred_spark_raw List containing per nucleotide Phred score values for each sequence
#' @slot phred_raw Numeric string containing mean Phred scores before trimming.
#' @slot phred_trim Numeric string containing mean Phred scores after trimming.
#' @slot Ns_raw Numeric string containing count of N's before trimming.
#' @slot Ns_trim Numeric string containing count of N's after trimming.
#' @slot length_raw Numeric string containing sequence length before trimming.
#' @slot length_trim Numeric string containing sequence length after trimming.
#' @slot seqs_raw Character string containing sequences before trimming.
#' @slot seqs_trim Character string containing sequence after trimming.
#' @slot decision Character string containing decision (PASS/FAIL) information based on \code{\link{isoQC}} 'min_phred_score' and 'min_length cutoffs'.
#' @slot input Character string containing input directory information.
#' @return Returns an class-isoQC object.
#' @seealso \code{\link{isoQC}}

setClass("isoQC",
         representation(
           date="character",
           trim.start.pos="numeric",
           trim.end.pos="numeric",
           filename="character",
           phred_spark_raw="list",
           phred_raw="numeric",
           phred_trim="numeric",
           Ns_raw="numeric",
           Ns_trim="numeric",
           length_raw="numeric",
           length_trim="numeric",
           seqs_raw="character",
           seqs_trim="character",
           decision="character",
           input="character"
         )
)

###########################################################################################################

#::::::::::::::::::::::::::::::::::::::::
# Set 'class-isoTAX' Class Object
#::::::::::::::::::::::::::::::::::::::::

#' @name class-isoTAX
#' @rdname class-isoTAX
#' @title isoTAX Class Object
#' @description S4 wrapper for \code{\link{isoTAX}} function. Access data via S4 slot functions.
#' @export
#' @slot input Character string containing input directory information.
#' @slot warning Character string containing list filenames of sequences that had poor alignment during taxonomic classification step.
#' @slot date Character string containing run date from each of the input Sanger sequence .ab1 files ("YYYY_MM_DD" format).
#' @slot filename Character string containing input filenames.
#' @slot phred_spark_raw List containing per nucleotide Phred score values for each sequence
#' @slot phred_raw Numeric string containing mean Phred scores before trimming.
#' @slot phred_trim Numeric string containing mean Phred scores after trimming.
#' @slot Ns_raw Numeric string containing count of N's before trimming.
#' @slot Ns_trim Numeric string containing count of N's after trimming.
#' @slot length_raw Numeric string containing sequence length before trimming.
#' @slot length_trim Numeric string containing sequence length after trimming.
#' @slot seqs_raw Character string containing sequences before trimming.
#' @slot seqs_trim Character string containing sequence after trimming.
#' @slot closest_match Character string containing species + type strain no. of closest match from reference database.
#' @slot NCBI_acc Character string containing NCBI accession number associated with closest match from reference database.
#' @slot ID Numeric string containing containing pairwise similarity value for query vs database reference sequence.
#' Calculation of ID is determined by isoTAX 'iddef' parameter (0-4, Default=2). See VSEARCH documentation for more details.
#'- (0) CD-HIT definition: (matching columns) / (shortest sequence length).
#'- (1) Edit distance: (matching columns) / (alignment length).
#'- (2) Edit distance excluding terminal gaps (default definition).
#'- (3) Marine Biological Lab definition counting each gap opening (internal or terminal) as a single mismatch, whether or not the gap was extended: 1.0- ((mismatches + gap openings)/(longest sequence length)).
#'- (4) BLAST definition, equivalent to --iddef 1 for global pairwise alignments.
#' @slot rank_phylum Character string containing Phylum rank taxonomy
#' @slot rank_class Character string containing Class rank taxonomy
#' @slot rank_order Character string containing Order rank taxonomy
#' @slot rank_family Character string containing Family rank taxonomy
#' @slot rank_genus Character string containing Genus rank taxonomy
#' @slot rank_species Character string containing Species rank taxonomy
#' @slot phylum_threshold Numeric string containing Phylum-level sequence similarity threshold for rank demarcation
#' @slot class_threshold Numeric string containing Class-level sequence similarity threshold for rank demarcation
#' @slot order_threshold Numeric string containing Order-level sequence similarity threshold for rank demarcation
#' @slot family_threshold Numeric string containing Family-level sequence similarity threshold for rank demarcation
#' @slot genus_threshold Numeric string containing Genus-level sequence similarity threshold for rank demarcation
#' @slot species_threshold Numeric string containing Species-level sequence similarity threshold for rank demarcation
#' @return Returns an class-isoTAX object.
#' @seealso \code{\link{isoTAX}}

#:::::::::::::::::::::::::::::::::::::
# Set class/method/generic functions
#:::::::::::::::::::::::::::::::::::::

setClass("isoTAX",
         representation(
           input="character",
           warning="character",
           date="character",
           filename="character",
           seqs_raw="character",
           phred_raw="numeric",
           Ns_raw="numeric",
           length_raw="numeric",
           phred_spark_raw="list",
           seqs_trim="character",
           phred_trim="numeric",
           Ns_trim="numeric",
           length_trim="numeric",
           closest_match="character",
           NCBI_acc="character",
           ID="numeric",
           rank_phylum="character",
           rank_class="character",
           rank_order="character",
           rank_family="character",
           rank_genus="character",
           rank_species="character",
           phylum_threshold="numeric",
           class_threshold="numeric",
           order_threshold="numeric",
           family_threshold="numeric",
           genus_threshold="numeric",
           species_threshold="numeric"
         )
)

###########################################################################################################

#::::::::::::::::::::::::::::::::::::::::
# Set 'class-isoLIB' Class Object
#::::::::::::::::::::::::::::::::::::::::

#' @name class-isoLIB
#' @rdname class-isoLIB
#' @title isoLIB Class Object
#' @description S4 wrapper for \code{\link{isoLIB}} function. Access data via S4 slot functions.
#' @export
#' @slot input Character string containing input directory information.
#' @slot sequence_group Character string containing list of group representative filenames.
#' @slot date Character string containing run date from each of the input Sanger sequence .ab1 files ("YYYY_MM_DD" format).
#' @slot filename Character string containing input filenames.
#' @slot phred_trim Numeric string containing mean Phred scores after trimming.
#' @slot Ns_trim Numeric string containing count of N's after trimming.
#' @slot length_trim Numeric string containing sequence length after trimming.
#' @slot seqs_trim Character string containing sequence after trimming.
#' @slot closest_match Character string containing species + type strain no. of closest match from reference database.
#' @slot NCBI_acc Character string containing NCBI accession number associated with closest match from reference database.
#' @slot ID Numeric string containing containing pairwise similarity value for query vs database reference sequence.
#' Calculation of ID is determined by isoTAX 'iddef' parameter (0-4, Default=2). See VSEARCH documentation for more details.
#'- (0) CD-HIT definition: (matching columns) / (shortest sequence length).
#'- (1) Edit distance: (matching columns) / (alignment length).
#'- (2) Edit distance excluding terminal gaps (default definition).
#'- (3) Marine Biological Lab definition counting each gap opening (internal or terminal) as a single mismatch, whether or not the gap was extended: 1.0- ((mismatches + gap openings)/(longest sequence length)).
#'- (4) BLAST definition, equivalent to --iddef 1 for global pairwise alignments.
#' @slot rank_phylum Character string containing Phylum rank taxonomy
#' @slot rank_class Character string containing Class rank taxonomy
#' @slot rank_order Character string containing Order rank taxonomy
#' @slot rank_family Character string containing Family rank taxonomy
#' @slot rank_genus Character string containing Genus rank taxonomy
#' @slot rank_species Character string containing Species rank taxonomy
#' @slot phylum_threshold Numeric string containing Phylum-level sequence similarity threshold for rank demarcation
#' @slot class_threshold Numeric string containing Class-level sequence similarity threshold for rank demarcation
#' @slot order_threshold Numeric string containing Order-level sequence similarity threshold for rank demarcation
#' @slot family_threshold Numeric string containing Family-level sequence similarity threshold for rank demarcation
#' @slot genus_threshold Numeric string containing Genus-level sequence similarity threshold for rank demarcation
#' @slot species_threshold Numeric string containing Species-level sequence similarity threshold for rank demarcation
#' @return Returns an class-isoLIB object.
#' @seealso \code{\link{isoLIB}}

setClass("isoLIB",
         representation(
           input="character",
           sequence_group="character",
           date="character",
           filename="character",
           seqs_trim="character",
           phred_trim="numeric",
           Ns_trim="numeric",
           length_trim="numeric",
           closest_match="character",
           NCBI_acc="character",
           ID="numeric",
           rank_phylum="character",
           rank_class="character",
           rank_order="character",
           rank_family="character",
           rank_genus="character",
           rank_species="character",
           representative="character",
           phylum_threshold="numeric",
           class_threshold="numeric",
           order_threshold="numeric",
           family_threshold="numeric",
           genus_threshold="numeric",
           species_threshold="numeric"
         )
)

###########################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::
# Set generic for 'isoQC' class object
#:::::::::::::::::::::::::::::::::::::::::::::::

#' @export
#' @aliases isoQC
setGeneric("isoQC",  function(input=NULL,
                              export_html=TRUE,
                              export_csv=TRUE,
                              export_fasta=TRUE,
                              export_fasta_revcomp=FALSE,
                              verbose=FALSE,
                              exclude=NULL,
                              min_phred_score=20,
                              min_length=200,
                              sliding_window_cutoff=NULL,
                              sliding_window_size = 15,
                              date=NULL,
                              files_manual=NULL) standardGeneric("isoQC"), signature=c("date"))


###########################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::
# Set generic for 'isoTAX' class object
#:::::::::::::::::::::::::::::::::::::::::::::::

#' @export
#' @aliases isoTAX
setGeneric("isoTAX",  function(input=NULL,
                               export_html=TRUE,
                               export_csv=TRUE,
                               export_blast_table=FALSE,
                               quick_search=FALSE,
                               db="16S_bac",
                               db_path=NULL,
                               vsearch_path=NULL,
                               iddef=2,
                               phylum_threshold = 75,
                               class_threshold = 78.5,
                               order_threshold = 82,
                               family_threshold = 86.5,
                               genus_threshold = 94.5,
                               species_threshold = 98.7) standardGeneric("isoTAX"), signature=c("input"))


###########################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::
# Set generic for 'isoLIB' class object
#:::::::::::::::::::::::::::::::::::::::::::::::

#' @export
#' @aliases isoLIB
setGeneric("isoLIB",  function(input=NULL,
                               old_lib_csv=NULL,
                               method="dark_mode",
                               group_cutoff = 0.995,
                               keep_old_reps=TRUE,
                               export_html=TRUE,
                               export_csv=TRUE,
                               include_warnings=TRUE,
                               vsearch_path=NULL,
                               phylum_threshold=75.0,
                               class_threshold=78.5,
                               order_threshold=82.0,
                               family_threshold=86.5,
                               genus_threshold=94.5,
                               species_threshold=98.7) standardGeneric("isoLIB"), signature=c("input"))


###########################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::
# Set generics for 'export_html' functions
#:::::::::::::::::::::::::::::::::::::::::::::::

#' @export
#' @aliases export_html
setGeneric("export_html", function(obj,
                                   method=NULL, 
                                   group_cutoff=NULL, 
                                   min_phred_score=NULL,
                                   min_length=NULL,
                                   sliding_window_cutoff=NULL,
                                   sliding_window_size = NULL,
                                   quick_search=NULL, db=NULL) standardGeneric("export_html"))




###########################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::
# Set generics for *_accessor functions
#:::::::::::::::::::::::::::::::::::::::::::::::


#' @export
setGeneric("trim.start.pos", function(obj) standardGeneric("trim.start.pos"))

#' @export
setGeneric("trim.start.pos<-", function(obj, value) standardGeneric("trim.start.pos<-"))




###########################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::
# Set method for 'isoQC' class object
#:::::::::::::::::::::::::::::::::::::::::::::::
#' @export
#' @name method-isoQC
#' @title setMethod functions for isoQC
#' @description Initiation of isoQC functions.
#' @rdname method-isoQC

setMethod("isoQC", signature(date="missing"), function(input=NULL,
                                                       export_html=TRUE,
                                                       export_csv=TRUE,
                                                       export_fasta=TRUE,
                                                       export_fasta_revcomp=FALSE,
                                                       verbose=FALSE,
                                                       exclude=NULL,
                                                       min_phred_score=20,
                                                       min_length=200,
                                                       sliding_window_cutoff=NULL,
                                                       sliding_window_size=15,
                                                       date=NULL,
                                                       files_manual=NULL) isoQC(input,
                                                                                export_html,
                                                                                export_csv,
                                                                                export_fasta,
                                                                                export_fasta_revcomp,
                                                                                verbose,
                                                                                exclude,
                                                                                min_phred_score,
                                                                                min_length,
                                                                                sliding_window_cutoff,
                                                                                sliding_window_size,
                                                                                date,
                                                                                files_manual))

#' @export
setMethod("trim.start.pos", "isoQC",
          function(obj) obj@trim.start.pos)

#' @export
setMethod("trim.start.pos<-", "isoQC",
          function(obj, value) {
            obj@trim.start.pos <- value
            obj
          })


###########################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::
# Set method for 'isoTAX' class object
#:::::::::::::::::::::::::::::::::::::::::::::::

#' @export
#' @name method-isoTAX
#' @title setMethod functions for isoTAX
#' @description Initiation of isoTAX functions.
#' @rdname method-isoTAX

setMethod("isoTAX", signature(input="missing"), function(input= NULL,
                                                         export_html=TRUE,
                                                         export_csv=TRUE,
                                                         export_blast_table=FALSE,
                                                         quick_search=FALSE,
                                                         db="16S_bac",
                                                         db_path=NULL,
                                                         vsearch_path=NULL,
                                                         iddef=2,
                                                         phylum_threshold = 75,
                                                         class_threshold = 78.5,
                                                         order_threshold = 82,
                                                         family_threshold = 86.5,
                                                         genus_threshold = 94.5,
                                                         species_threshold = 98.7) isoTAX(input,
                                                                                       export_html,
                                                                                       export_csv,
                                                                                       export_blast_table,
                                                                                       quick_search,
                                                                                       db,
                                                                                       db_path,
                                                                                       vsearch_path,
                                                                                       iddef,
                                                                                       phylum_threshold,
                                                                                       class_threshold ,
                                                                                       order_threshold,
                                                                                       family_threshold,
                                                                                       genus_threshold,
                                                                                       species_threshold))


###########################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::
# Set method for 'isoLIB' class object
#:::::::::::::::::::::::::::::::::::::::::::::::

#' @export
#' @name method-isoLIB
#' @title setMethod functions for isoLIB
#' @description Initiation of isoLIB functions.
#' @rdname method-isoLIB
setMethod("isoLIB", signature(input="missing"), function(input=NULL,
                                                         old_lib_csv=NULL,
                                                         method="dark_mode",
                                                         group_cutoff = 0.995,
                                                         keep_old_reps=TRUE,
                                                         export_html=TRUE,
                                                         export_csv=TRUE,
                                                         include_warnings=TRUE,
                                                         vsearch_path=NULL,
                                                         phylum_threshold=75.0,
                                                         class_threshold=78.5,
                                                         order_threshold=82.0,
                                                         family_threshold=86.5,
                                                         genus_threshold=94.5,
                                                         species_threshold=98.7) isoLIB(input,
                                                                                     old_lib_csv,
                                                                                     method,
                                                                                     group_cutoff,
                                                                                     keep_old_reps,
                                                                                     export_html,
                                                                                     export_csv,
                                                                                     include_warnings,
                                                                                     vsearch_path,
                                                                                     phylum_threshold,
                                                                                     class_threshold,
                                                                                     order_threshold,
                                                                                     family_threshold,
                                                                                     genus_threshold,
                                                                                     species_threshold))


###########################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::
# Set methods for 'show' S4 generic
#:::::::::::::::::::::::::::::::::::::::::::::::

#-------------------------------
#' @export
#' @name show
#' @title Generic show method for S4 class objects
#' @description Generic show method for S4 class objects.
#' @rdname show
#' @import methods

setMethod("show", "isoQC",
          function(object) {
            cat(paste("isolateR library file overview:"), "\n")
            cat(paste("----------------------------------------------------------------------------------"), "\n")
            cat(paste("Total sequences = ", length(object@date), sep=""), "\n")
            cat(paste("Mean Phred score before/after trimming = \t|   ", round(mean(object@phred_raw),2),"\t |   ",round(mean(object@phred_trim),2),"\t |", sep=""), "\n")
            cat(paste("Mean sequence length before/after trimming = \t|   ", round(mean(object@length_raw),2),"\t |   ",round(mean(object@length_trim),2), "\t |", sep=""), "\n")
            cat(paste("Mean number of N's before/after trimming = \t|   ", format(round(mean(object@Ns_raw),2),nsmall=2) ,"\t |   ",format(round(mean(object@Ns_trim),2),nsmall=2) , "\t |", sep=""), "\n")
          }
)

#-------------------------------
#' @export
#' @name show
#' @rdname show
#' @import methods

setMethod("show", "isoTAX",
          function(object) {
            cat(paste("isolateR library file overview:"), "\n")
            cat(paste("========================================================="), "\n")
            cat(paste("Total sequences                             = \t   ", length(object@date), sep=""), "\n")
            cat(paste("Number of sequences under Phylum threshod   = \t   ", length(isoTAX.S4@ID[isoTAX.S4@ID < isoTAX.S4@phylum_threshold[1]]), sep=""), "\n")
            cat(paste("Number of sequences under Class threshod    = \t   ", length(isoTAX.S4@ID[isoTAX.S4@ID < isoTAX.S4@class_threshold[1]]), sep=""), "\n")
            cat(paste("Number of sequences under Order threshod    = \t   ", length(isoTAX.S4@ID[isoTAX.S4@ID < isoTAX.S4@order_threshold[1]]), sep=""), "\n")
            cat(paste("Number of sequences under Family threshod   = \t   ", length(isoTAX.S4@ID[isoTAX.S4@ID < isoTAX.S4@family_threshold[1]]), sep=""), "\n")
            cat(paste("Number of sequences under Genus threshod    = \t   ", length(isoTAX.S4@ID[isoTAX.S4@ID < isoTAX.S4@genus_threshold[1]]), sep=""), "\n")
            cat(paste("Number of sequences under Species threshod  = \t   ", length(isoTAX.S4@ID[isoTAX.S4@ID < isoTAX.S4@species_threshold[1]]), sep=""), "\n")
          }
)
#-------------------------------


#' @export
#' @name show
#' @rdname show
#' @import methods

setMethod("show", "isoLIB",
          function(object) {
            cat(paste("isolateR - isoLIB class file overview:"), "\n")
            cat(paste("----------------------------------------------------------------------------------"), "\n")
            cat(paste("Total sequences = ", length(object@date), sep=""), "\n")
          }
)



###########################################################################################################
#:::::::::::::::::::::::::::::::::::::::::::::::
# Set S4 conversion functions
#:::::::::::::::::::::::::::::::::::::::::::::::

#-------------------------------
#' @export
#' @rdname S4_to_dataframe
#' @title Converts S4 objects (isoQC, isoTAX, or isoLIB) to dataframe
#' @description Helper function to convert S4 class objects (\code{\link{isoQC}}, \code{\link{isoTAX}}, or \code{\link{isoLIB}}) to dataframe
#' @param obj S4 object generated from \code{\link{isoQC}}, \code{\link{isoTAX}}, or \code{\link{isoLIB}} steps
#' @returns Returns a dataframe containing sequence information in columns.
#' @import dplyr
#' @import methods

S4_to_dataframe <- function(obj) {
  nms <- slotNames(obj)
  lst <- sapply(nms, function(nm) slot(obj, nm))
  if("phred_spark_raw" %in% nms ==TRUE){
    lst$phred_spark_raw <- sapply(lst$phred_spark_raw, function(x) paste(unlist(x), collapse="_"))
  }
  df.out <- as.data.frame(lst)
  if(length(obj@filename) < 2){
    df.out <- df.out[1,]
    if(class(obj)[1]!="isoLIB"){
      df.out$phred_spark_raw <- paste(unlist(lst$phred_spark_raw), collapse="_")
    }
    if(class(obj)[1]=="isoLIB"){
      df.out <- as.data.frame(t(lst))
    }
  }
  return(df.out)
}

#-------------------------------
#' @export
#' @rdname df_to_isoTAX
#' @title Convert isoTAX .CSV output to isoTAX class object
#' @description Helper function to convert isoTAX .CSV output to a \code{\link{class-isoTAX}} class object.
#' @param df Dataframe in same format as .CSV output file from \code{\link{isoTAX}} step.
#' @returns Returns an S4 \code{\link{class-isoTAX}} object that can be used to generate interactive HTML output tables.
#' @import dplyr
#' @import methods


df_to_isoTAX <- function(df) {
  phred_spark_raw <- NULL
  . <- NULL
  isotax.S4 <- new("isoTAX", input=getwd())
  # Make phred scores into lists for sparkline formatting (V1=raw, V2=trim)
  isotax.S4@warning <- df$warning
  isotax.S4@date <- df$date
  isotax.S4@filename <- df$filename
  isotax.S4@seqs_raw <- df$seqs_raw
  isotax.S4@phred_raw <- df$phred_raw
  isotax.S4@Ns_raw <- df$Ns_raw
  isotax.S4@length_raw <- df$length_raw
  isotax.S4@phred_spark_raw <- (df %>% mutate(phred_spark_raw = strsplit(phred_spark_raw, '_')) %>%  mutate(phred_spark_raw = lapply(.$phred_spark_raw, function(x) list(as.numeric(x)))))$phred_spark_raw
  isotax.S4@seqs_trim <- df$seqs_trim
  isotax.S4@phred_trim <- df$phred_trim
  isotax.S4@Ns_trim <- df$Ns_trim
  isotax.S4@length_trim <- df$length_trim
  isotax.S4@closest_match <- df$closest_match
  isotax.S4@NCBI_acc <- df$NCBI_acc
  isotax.S4@ID <- df$ID
  isotax.S4@rank_phylum <- df$rank_phylum
  isotax.S4@rank_class <- df$rank_class
  isotax.S4@rank_order <- df$rank_order
  isotax.S4@rank_family <- df$rank_family
  isotax.S4@rank_genus <- df$rank_genus
  isotax.S4@rank_species <- df$rank_species
  
  return(isotax.S4)
  
}

#-------------------------------
#' @export
#' @rdname df_to_isoLIB
#' @title Convert isoLIB .CSV output to isoLIB class object
#' @description Helper function to convert isoLIB .CSV output to a \code{\link{class-isoLIB}} class object.
#' @param df Dataframe in same format as .CSV output file from \code{\link{isoLIB}} step.
#' @returns Returns an S4 \code{\link{class-isoLIB}} object that can be used to generate interactive HTML output tables.
#' @import methods


df_to_isoLIB <- function(df) {
  
  isolib.S4 <- new("isoLIB", input=getwd())
  # Make phred scores into lists for sparkline formatting (V1=raw, V2=trim)
  isolib.S4@sequence_group <- df$sequence_group
  isolib.S4@filename <- df$filename
  isolib.S4@date <- df$date
  isolib.S4@seqs_trim <- df$seqs_trim
  isolib.S4@phred_trim <- df$phred_trim
  isolib.S4@Ns_trim <- df$Ns_trim
  isolib.S4@length_trim <- df$length_trim
  isolib.S4@closest_match <- df$closest_match
  isolib.S4@NCBI_acc <- df$NCBI_acc
  isolib.S4@ID <- df$ID
  isolib.S4@rank_phylum <- df$rank_phylum
  isolib.S4@rank_class <- df$rank_class
  isolib.S4@rank_order <- df$rank_order
  isolib.S4@rank_family <- df$rank_family
  isolib.S4@rank_genus <- df$rank_genus
  isolib.S4@rank_species <- df$rank_species
  isolib.S4@representative <- df$representative
  isolib.S4@phylum_threshold <- df$phylum_threshold
  isolib.S4@class_threshold <- df$class_threshold
  isolib.S4@order_threshold <- df$order_threshold
  isolib.S4@family_threshold <- df$family_threshold
  isolib.S4@genus_threshold <- df$genus_threshold
  isolib.S4@species_threshold <- df$species_threshold
  return(isolib.S4)
}

###########################################################################################################

#:::::::::::::::::::::::::::::::::::::::::::::::
# Set methods for 'export_html' S4 generic
#:::::::::::::::::::::::::::::::::::::::::::::::
#' @export
#' @name export_html
#' @rdname export_html
#' @aliases export_html-isoQC
#' @title Export HTML for isoQC > isoTAX > isoLIB class objects
#' @description S4 wrapper functions to export interactive HTML tables from \code{\link{isoQC}}, \code{\link{isoTAX}}, or \code{\link{isoLIB}} class objects. Saves to HTML to current working directory and automatically opens.
#' @param obj An S4 class object generated from one of \code{\link{isoQC}}, \code{\link{isoTAX}}, or \code{\link{isoLIB}} steps
#' @returns HTML output file saved to working directory.
#' @importFrom crosstalk bscols
#' @importFrom crosstalk filter_checkbox
#' @importFrom crosstalk filter_select
#' @importFrom crosstalk filter_slider
#' @importFrom crosstalk SharedData
#' @importFrom pander openFileInOS
#' @importFrom dataui dui_sparkline
#' @importFrom dataui dui_sparkbandline
#' @importFrom dataui dui_sparklineseries
#' @importFrom htmltools save_html
#' @importFrom htmltools div
#' @importFrom htmltools browsable
#' @import reactable
#' @import reactablefmtr
#' @importFrom scales comma
#' @importFrom scales comma_format
#' @importFrom scales label_number

#:::::::::::::::::::::::::::::::
# MAKE HTML #1: QC TABLE
#:::::::::::::::::::::::::::::::

setMethod("export_html", "isoQC",
          function(obj, min_phred_score=NULL, min_length=NULL, sliding_window_cutoff=NULL, sliding_window_size = NULL){
            #Import checks
            if(is.null(min_phred_score)){ message("The 'min_phred_score' parameter is blank. Setting to 1 to enable graphing.")
              min_phred_score <- 1}
            if(is.null(min_length)){ message("The 'min_length' parameter is blank. Setting to 1 to enable graphing.")
              min_length <- 1}
            if(is.null(sliding_window_cutoff)) {sliding_window_cutoff <- "Auto" }
            if(is.null(sliding_window_size)) { message("The 'sliding_window_size' parameter is blank. Setting to 1 to enable graphing.")
              sliding_window_size <- 1}
            
            #Import S4 object
            isoQC.df <- S4_to_dataframe(obj) %>%
              mutate(phred_spark_raw = strsplit(phred_spark_raw, '_')) %>%
              mutate(phred_spark_raw = lapply(.$phred_spark_raw, function(x) list(as.numeric(x)))) %>%
              select(date,
                     filename,
                     seqs_raw,
                     phred_raw,
                     Ns_raw,
                     length_raw,
                     phred_spark_raw,
                     seqs_trim,
                     phred_trim,
                     Ns_trim,
                     length_trim,
                     decision)
            
            suppressWarnings(isoQC.df.g <- rbind(isoQC.df %>% select(filename, phred_raw, Ns_raw, length_raw, decision) %>% mutate(group="raw") %>% `colnames<-`(c("filename", "phred", "Ns", "length", "decision", "group")),
                                                 isoQC.df %>% select(filename, phred_trim, Ns_trim, length_trim, decision) %>% mutate(group="trim") %>% `colnames<-`(c("filename", "phred", "Ns", "length", "decision", "group"))) %>%
                               mutate_at(vars(decision), funs(factor(., levels=c("Pass", "Fail")))))
            
            #-----------------------------------------------------------------------------------------------------------Summary plots
            g.phred <- ggplot2::ggplot(isoQC.df.g, ggplot2::aes(x=group, y=phred, colour=decision)) +
              ggplot2::geom_hline(yintercept = min_phred_score, linewidth=1.5, colour="red", alpha=0.2, linetype="solid") +
              ggplot2::annotate(geom="text", vjust=1.75, size=3.3, x=0.7, y=min_phred_score, label=paste('min_phred_score="', min_phred_score, '"', sep=""), color=scales::alpha("red", 0.85)) +  
              ggiraph::geom_boxplot_interactive(ggplot2::aes(fill=decision, color=decision, #data_id = group,
                                                             tooltip = ggplot2::after_stat({
                                                               paste0(
                                                                 "Decision: ", .data$fill,
                                                                 "\nQ1: ", prettyNum(.data$ymin),
                                                                 "\nQ3: ", prettyNum(.data$ymax),
                                                                 "\nmedian: ", prettyNum(.data$middle)
                                                               )
                                                             })), position = ggplot2::position_dodge(width = 0.8), outlier.shape = NA) +
              ggiraph::geom_point_interactive(ggplot2::aes(group=decision, colour=decision, data_id = filename, tooltip=filename), position=ggplot2::position_jitterdodge(dodge.width = 0.8), size=0.7, alpha=0.5, hover_nearest = TRUE) +
              ggplot2::scale_color_manual(values=c(Pass=scales::alpha("black", 0.7),Fail=scales::alpha("black", 0.7)), name="Decision") + 
              ggplot2::scale_fill_manual(values=c(Pass=scales::alpha("#CCE6CC", 0.9), Fail=scales::alpha("#B20A2C", 0.5)  ), name="Decision") +
              ggplot2::geom_vline(xintercept = c(1.5), linewidth=1, colour="black", linetype="dotted") +
              ggplot2::scale_x_discrete(labels=c("raw"="Raw","trim"="Trimmed"))  +
              ggplot2::ylab("Phred Quality Score") + 
              ggplot2::theme(panel.background = ggplot2::element_blank(),
                             panel.border = ggplot2::element_rect(color ="black", fill = NA, linewidth = 0.5, linetype = 1),
                             panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                             legend.title = ggplot2::element_text(size=16),legend.text = ggplot2::element_text(size=12),
                             strip.background = ggplot2::element_blank(), strip.text.x = ggplot2::element_blank(),
                             axis.title.x = ggplot2::element_blank(), axis.text.x=ggplot2::element_text(size=12), axis.ticks.x=ggplot2::element_blank(),
                             axis.title.y = ggplot2::element_text(size=16), axis.text.y=ggplot2::element_text(size=12))
            
            g.length <- ggplot2::ggplot(isoQC.df.g, ggplot2::aes(x=group, y=length, colour=decision)) + 
              ggplot2::geom_hline(yintercept = min_length, linewidth=1.5, colour="red", alpha=0.2, linetype="solid") +
              ggplot2::annotate(geom="text", vjust=1.75, size=3.3, x=0.65, y=min_length, label=paste('min_length="', min_length, '"', sep=""), color=scales::alpha("red", 0.85)) +
              ggiraph::geom_boxplot_interactive(ggplot2::aes(fill=decision, #data_id = group,
                                                             tooltip = ggplot2::after_stat({
                                                               paste0(
                                                                 "Decision: ", .data$fill,
                                                                 "\nQ1: ", prettyNum(.data$ymin),
                                                                 "\nQ3: ", prettyNum(.data$ymax),
                                                                 "\nmedian: ", prettyNum(.data$middle)
                                                               )
                                                             })), position = ggplot2::position_dodge(width = 0.8), outlier.shape = NA) +
              ggiraph::geom_point_interactive(ggplot2::aes(group=decision, colour=decision, data_id = filename, tooltip=filename), position=ggplot2::position_jitterdodge(dodge.width = 0.8), size = 0.7, alpha=0.5, hover_nearest = TRUE) +
              ggplot2::scale_color_manual(values=c(Pass=scales::alpha("black", 0.7),Fail=scales::alpha("black", 0.7)),  name="Decision") + 
              ggplot2::scale_fill_manual(values=c(Pass=scales::alpha("#CCE6CC", 0.9), Fail=scales::alpha("#B20A2C", 0.5)  ), name="Decision") +
              ggplot2::geom_vline(xintercept = c(1.5), linewidth=1, colour="black", linetype="dotted") +
              ggplot2::scale_x_discrete(labels=c("raw"="Raw","trim"="Trimmed"))  +
              ggplot2::ylab("Sequence length (nt)") + 
              ggplot2::theme(panel.background = ggplot2::element_blank(),
                             panel.border = ggplot2::element_rect(color ="black", fill = NA, linewidth = 0.5, linetype = 1),
                             panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                             legend.title = ggplot2::element_text(size=16),legend.text = ggplot2::element_text(size=12),
                             strip.background = ggplot2::element_blank(), strip.text.x = ggplot2::element_blank(),
                             axis.title.x = ggplot2::element_blank(), axis.text.x=ggplot2::element_text(size=12), axis.ticks.x=ggplot2::element_blank(),
                             axis.title.y = ggplot2::element_text(size=16), axis.text.y=ggplot2::element_text(size=12))
            
            isoQC.plots1 <- cowplot::plot_grid(g.phred + ggplot2::theme(legend.position="none"),
                                               NULL,
                                               g.length + ggplot2::theme(legend.position="none"),
                                               suppressWarnings(cowplot::get_legend(g.phred)), 
                                               ncol=4, rel_widths=c(1,0.1,1,0.25))
            
            ggiraph::set_girafe_defaults(
              opts_hover_inv = ggiraph::opts_hover_inv(css = "opacity:0.4;"), 
              opts_hover = ggiraph::opts_hover(css = "fill:orange;stroke:orange;stroke-width:7;"),
              #opts_zoom = ggiraph::opts_zoom(min = 1, max = 4),
              opts_tooltip = ggiraph::opts_tooltip(css = "padding:3px;background-color:#333333;color:white;"),
              opts_sizing = ggiraph::opts_sizing(rescale = FALSE),
              opts_toolbar = ggiraph::opts_toolbar(saveaspng = TRUE, pngname="isoQC", delay_mouseout = 5000, position="topright", hidden=c("selection", "zoom"))
            )
            
            isoQC.plots <- ggiraph::girafe(ggobj = isoQC.plots1, width_svg = 14.1, height_svg = 4.5)
            
            
            #-----------------------------------------------------------------------------------------------------------
            
            #Set universal font size
            uni.font <- 11
            
            suppressWarnings(
              html_output <- crosstalk::bscols(
                reactable::reactable(isoQC.df,
                                     fullWidth=FALSE,
                                     searchable = TRUE,
                                     bordered = TRUE,
                                     resizable =TRUE,
                                     defaultPageSize=20,
                                     pageSizeOptions = c(10, 25, 50, 100, 1000),
                                     showPageSizeOptions = TRUE,
                                     highlight = TRUE,
                                     showSortable = TRUE,
                                     compact=TRUE,
                                     wrap = TRUE,
                                     theme = reactableTheme(
                                       headerStyle = list(
                                         fontFamily = "sans-serif", fontWeight="bold", fontSize=uni.font,
                                         "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                                         "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)")),
                                       style = list(
                                         fontFamily = "sans-serif", fontWeight="normal",
                                         "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                                         "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"))),
                                     style = list(minWidth = 1400),
                                     #--------------------Universal settings
                                     defaultColDef = colDef(minWidth= 75,
                                                            align="center",
                                                            style = list(whiteSpace = "nowrap",
                                                                         fontSize=uni.font,
                                                                         align="center",
                                                                         fontFamily = "sans-serif")),
                                     #--------------------
                                     columns = list(
                                       #--------------------
                                       # Spark lines
                                       #-------------------
                                       phred_spark_raw = colDef(name = "Quality sparkline", width = 200, cell = function(value, index) {
                                         dataui::dui_sparkline(data = value[[1]],
                                                               height = 25,
                                                               min = min(unlist(isoQC.df$phred_spark_raw)),
                                                               max = max(unlist(isoQC.df$phred_spark_raw)),
                                                               margin = list( top= 2,
                                                                              right= 5,
                                                                              bottom= 2,
                                                                              left= 5 ),
                                                               components = list(
                                                                 dataui::dui_sparkbandline(band = list( from = list( x = obj@trim.start.pos[[index]] ),
                                                                                                        to = list( x = obj@trim.end.pos[[index]] ) ),
                                                                                           fill = "green",
                                                                                           fillOpacity=0.2 ), #, fillOpacity=0.3 #use if fill not banding "url(#band_pattern_misc)"
                                                                 dataui::dui_sparklineseries(strokeWidth = 0.5,
                                                                                             stroke = "black")
                                                               )
                                         )
                                       }
                                       ),
                                       date = colDef(name = "Date", width=75, style=list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                       seqs_raw = colDef(name = "rawSeq",
                                                         style=list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                       seqs_trim = colDef(name = "trimSeq",
                                                          style=list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                       filename = colDef(minWidth= 125,
                                                         name ="Filename",
                                                         style = list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                       decision = colDef(name="Decision", vAlign = "center",  minWidth=170,
                                                         style = function(value, index, decision){
                                                           if(isoQC.df$decision[[index]] == "Pass") {
                                                             color.background <- grDevices::adjustcolor("#CCE6CC", alpha.f=1)
                                                             color.text <- "black"
                                                           } else if(!isoQC.df$decision[[index]] == "Pass") {
                                                             color.background <- grDevices::adjustcolor("#B20A2C", alpha.f=0.7)
                                                             color.text <- "white" }
                                                           list(background = color.background, color = color.text,
                                                                whiteSpace = "nowrap",
                                                                fontSize=uni.font+6,
                                                                align="center",
                                                                fontFamily = "sans-serif")}),
                                       decision_col = colDef(show=FALSE),
                                       phred_raw = colDef(name = "rawQ",
                                                          minWidth= 75,
                                                          cell = gauge_chart(isoQC.df,
                                                                             fill_color = c('#e5f5e0', '#a1d99b', '#31a354'),
                                                                             opacity = 0.5,
                                                                             number_fmt = scales::comma,
                                                                             text_size = 12,
                                                                             max_value = max(c(isoQC.df$phred_raw,isoQC.df$phred_trim)),
                                                                             show_min_max = FALSE)),
                                       phred_trim = colDef(name = "trimQ",
                                                           minWidth= 75,
                                                           cell = gauge_chart(isoQC.df,
                                                                              fill_color = c('#e5f5e0', '#a1d99b', '#31a354'),
                                                                              opacity = 0.5,
                                                                              number_fmt = scales::comma,
                                                                              #bold_text = TRUE,
                                                                              text_size = 12,
                                                                              max_value = max(c(isoQC.df$phred_raw,isoQC.df$phred_trim)),
                                                                              show_min_max = FALSE)),
                                       Ns_raw = colDef(name = "rawNs",
                                                       minWidth= 75,
                                                       cell = gauge_chart(isoQC.df,
                                                                          fill_color = c("white", "#B20A2C"),
                                                                          opacity = 0.5,
                                                                          number_fmt = scales::comma,
                                                                          text_size = 12,
                                                                          max_value = max(c(isoQC.df$Ns_raw,isoQC.df$Ns_trim)),
                                                                          show_min_max = FALSE)),
                                       Ns_trim = colDef(name = "trimNs",
                                                        minWidth= 75,
                                                        cell = gauge_chart(isoQC.df,
                                                                           fill_color = c("white", "#B20A2C"),
                                                                           opacity = 0.5,
                                                                           number_fmt = scales::comma,
                                                                           text_size = 12,
                                                                           max_value = max(c(isoQC.df$Ns_raw,isoQC.df$Ns_trim)),
                                                                           show_min_max = FALSE)),
                                       # Barplot widgets
                                       #-------------------
                                       length_raw = colDef(name = "rawLength",
                                                           minWidth=150,
                                                           cell = data_bars(isoQC.df,
                                                                            text_position = "inside-base",
                                                                            text_size = 12,
                                                                            fill_color = c('#FFF2D9','#FFE1A6','#FFCB66','#FFB627'),
                                                                            fill_gradient = TRUE,
                                                                            background = 'transparent',
                                                                            number_fmt = scales::comma_format(accuracy = 0.1),
                                                                            round_edges = FALSE, align_bars="left")),
                                       length_trim = colDef(name = "trimLength",
                                                            minWidth=150,
                                                            cell = data_bars(isoQC.df,
                                                                             text_position = "inside-base", #inside-end / outside-base / above
                                                                             text_size = 12,
                                                                             #fill_color="#6f86ab",
                                                                             fill_color = c('#FFF2D9','#FFE1A6','#FFCB66','#FFB627'),
                                                                             fill_gradient = TRUE,
                                                                             background = 'transparent',
                                                                             number_fmt = scales::comma_format(accuracy = 0.1),
                                                                             round_edges = FALSE, align_bars="left")) )) %>% reactablefmtr::google_font(font_family = "Source Sans Pro") %>%
                  htmlwidgets::prependContent(
                    shiny::tags$div(
                      shiny::tags$div(
                        shiny::tags$div(
                          "isoQC output table", 
                          style = htmltools::css('font-size' = '24pt', 'font-weight' = 'bold',
                                                 'text-align' = 'left', 'margin-bottom' = 0, 'padding-top' = '10px', 'padding-left' = '10px', 'vertical-align' = 'middle')),
                        shiny::tags$div(paste("------------------------------------------------------------", sep=""), 
                                        style = htmltools::css('font-size' = '14pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste('Project directory: "', unlist(strsplit(obj@input, '/'))[length(unlist(strsplit(obj@input, '/')))], '"', sep=""), 
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'bold', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("Date last updated: ", Sys.Date(), sep=""), 
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("------------------------------------------------------------", sep=""), 
                                        style = htmltools::css('font-size' = '14pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("Settings:", sep=""), 
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'bold', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("min_phred_score = '", min_phred_score,"'", sep=""),
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-top' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("min_length = '", min_length,"'", sep=""),
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("sliding_window_cutoff = '", sliding_window_cutoff,"'", sep=""), 
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("sliding_window_size = '", sliding_window_size,"'",sep=""), 
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("Results:", sep=""), 
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'bold', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '10px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("Total sequences:       ", nrow(isoQC.df), sep=""), 
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("Mean Phred quality before/after trimming:     ", round(mean(isoQC.df$phred_raw)), "\t|  ", round(mean(isoQC.df$phred_trim)), sep=""),
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("Mean number of Ns before/after trimming:      ", round(mean(isoQC.df$Ns_raw)), "\t|  ", round(mean(isoQC.df$Ns_trim)), sep=""),
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        shiny::tags$div(paste("Mean sequence length before/after trimming:   ", round(mean(isoQC.df$length_raw)), "\t|  ", round(mean(isoQC.df$length_trim)), sep=""),
                                        style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                               'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                        style = htmltools::css(width = '28%')),
                      shiny::tags$div(isoQC.plots,
                                      style = htmltools::css(width = '72%')),
                      style = htmltools::css(
                        width = '1400px',
                        display = 'inline-flex')))
              )
            )
            
            path <- obj@input
            output <- "isolateR_output/01_isoQC_results.html"
            fname_html <- file.path(path, output)
            suppressMessages(htmltools::save_html(html_output, fname_html))
            writeLines(gsub('<meta charset="utf-8"/>', paste('<meta charset="utf-8"/>\n<title>isoQC > ', unlist(strsplit(obj@input, '/'))[length(unlist(strsplit(obj@input, '/')))],'</title>', sep=""), readLines(fname_html)), fname_html)
            pander::openFileInOS(fname_html)
            #return(html_output)
          })

#:::::::::::::::::::::::::::::::
# MAKE HTML #2: TAX TABLE
#:::::::::::::::::::::::::::::::
#' @export
#' @name export_html
#' @rdname export_html
#' @aliases export_html-isoTAX
#' @importFrom crosstalk bscols
#' @importFrom crosstalk filter_checkbox
#' @importFrom crosstalk filter_select
#' @importFrom crosstalk filter_slider
#' @importFrom crosstalk SharedData
#' @importFrom pander openFileInOS
#' @importFrom dataui dui_sparkline
#' @importFrom dataui dui_sparkbandline
#' @importFrom dataui dui_sparklineseries
#' @importFrom htmltools save_html
#' @importFrom htmltools div
#' @importFrom htmltools browsable
#' @import reactable
#' @import reactablefmtr
#' @importFrom scales comma
#' @importFrom scales comma_format
#' @importFrom scales label_number


setMethod("export_html", "isoTAX",
          function(obj, quick_search=NULL, db=NULL){
            #Import checks
            if(is.null(quick_search)){ message("The 'quick_search' parameter is blank. Setting to '' to enable graphing.")
              quick_search <- ''}
            if(is.null(db)){ message("The 'db' parameter is blank. Setting to '' to enable graphing.")
              db <- ''}
            
            #Import S4 object
            merged_input <- S4_to_dataframe(obj) %>%
              mutate(phred_spark_raw = strsplit(phred_spark_raw, '_')) %>%
              mutate(phred_spark_raw = lapply(.$phred_spark_raw, function(x) list(as.numeric(x)))) %>%
              select(-input, -phylum_threshold, -class_threshold, -order_threshold, -family_threshold, -genus_threshold, -species_threshold)
            path <- getwd()
            output="02_isoTAX_results.html"
            
            #Set taxonomic rank cutoffs
            phylum_threshold <- obj@phylum_threshold
            class_threshold <- obj@class_threshold
            order_threshold <- obj@order_threshold
            family_threshold <- obj@family_threshold
            genus_threshold <- obj@genus_threshold
            species_threshold <- obj@species_threshold
            #-----------------------------------------------------------------------------------------------------------------Make sunplot
            
            d1 <- data.frame(labels = gsub(".", "_", merged_input$filename, fixed=TRUE),
                             parents = merged_input$rank_species) %>% mutate(ids = paste(labels, sep=""))
            d2 <- data.frame(labels = stringr::str_split_fixed(unique(paste(merged_input$rank_species, merged_input$rank_genus, sep="___")), "___", 2)[,1],
                             parents = stringr::str_split_fixed(unique(paste(merged_input$rank_species, merged_input$rank_genus, sep="___")), "___", 2)[,2]) %>% mutate(ids = paste(labels, sep=""))
            d3 <- data.frame(labels = stringr::str_split_fixed(unique(paste(merged_input$rank_genus, merged_input$rank_family, sep="___")), "___", 2)[,1],
                             parents = stringr::str_split_fixed(unique(paste(merged_input$rank_genus, merged_input$rank_family, sep="___")), "___", 2)[,2]) %>% mutate(ids = paste(labels, sep=""))
            d4 <- data.frame(labels = stringr::str_split_fixed(unique(paste(merged_input$rank_family, merged_input$rank_order, sep="___")), "___", 2)[,1],
                             parents = stringr::str_split_fixed(unique(paste(merged_input$rank_family, merged_input$rank_order, sep="___")), "___", 2)[,2]) %>% mutate(ids = paste(labels, sep=""))
            d5 <- data.frame(labels = stringr::str_split_fixed(unique(paste(merged_input$rank_order, merged_input$rank_class, sep="___")), "___", 2)[,1],
                             parents = stringr::str_split_fixed(unique(paste(merged_input$rank_order, merged_input$rank_class, sep="___")), "___", 2)[,2]) %>% mutate(ids = paste(labels, sep=""))
            d6 <- data.frame(labels = stringr::str_split_fixed(unique(paste(merged_input$rank_class, merged_input$rank_phylum, sep="___")), "___", 2)[,1],
                             parents = stringr::str_split_fixed(unique(paste(merged_input$rank_class, merged_input$rank_phylum, sep="___")), "___", 2)[,2]) %>% mutate(ids = paste(labels, sep=""))
            
            d2e <- data.frame(labels = unique(gsub(" ", "_", merged_input$rank_phylum)),parents = paste("Phlyum", sep="")) %>% mutate(ids = paste(labels, sep=""))
            d1e <- data.frame(labels = c("Phlyum"), parents = c("")) %>% mutate(ids = paste(labels, sep=""))
            sunplot <- rbind(d1e, d2e, d6, d5, d4, d3, d2, d1) %>% select(ids, labels, parents) %>% mutate(labels = gsub(" ", "<br>", labels))
            
            if(n_distinct(merged_input$rank_class) <10 ){
              d2e <- data.frame(labels = unique(gsub(" ", "_", merged_input$rank_class)),parents = paste("Class", sep="")) %>% mutate(ids = paste(labels, sep=""))
              d1e <- data.frame(labels = c("Class"), parents = c("")) %>% mutate(ids = paste(labels, sep=""))
              sunplot <- rbind(d1e, d2e, d5, d4, d3, d2, d1) %>% select(ids, labels, parents) %>% mutate(labels = gsub(" ", "<br>", labels))
            }
            
            if(n_distinct(merged_input$rank_order) <10 ){
              d2e <- data.frame(labels = unique(gsub(" ", "_", merged_input$rank_order)),parents = paste("Order", sep="")) %>% mutate(ids = paste(labels, sep=""))
              d1e <- data.frame(labels = c("Order"), parents = c("")) %>% mutate(ids = paste(labels, sep=""))
              sunplot <- rbind(d1e, d2e, d4, d3, d2, d1) %>% select(ids, labels, parents) %>% mutate(labels = gsub(" ", "<br>", labels))
            }
            if(n_distinct(merged_input$rank_family) <10 ){
              d2e <- data.frame(labels = unique(gsub(" ", "_", merged_input$rank_family)),parents = paste("Family", sep="")) %>% mutate(ids = paste(labels, sep=""))
              d1e <- data.frame(labels = c("Family"), parents = c("")) %>% mutate(ids = paste(labels, sep=""))
              sunplot <- rbind(d1e, d2e, d3, d2, d1) %>% select(ids, labels, parents) %>% mutate(labels = gsub(" ", "<br>", labels))
            }
            if(n_distinct(merged_input$rank_genus) <10 ){
              d2e <- data.frame(labels = unique(gsub(" ", "_", merged_input$rank_genus)),parents = paste("Genus", sep="")) %>% mutate(ids = paste(labels, sep=""))
              d1e <- data.frame(labels = c("Genus"), parents = c("")) %>% mutate(ids = paste(labels, sep=""))
              sunplot <- rbind(d1e, d2e, d2, d1) %>% select(ids, labels, parents) %>% mutate(labels = gsub(" ", "<br>", labels))
            }
            
            sunplot_small <- plotly::plot_ly(sunplot, ids = ~ids, labels = ~labels, parents = ~parents, type = 'sunburst', width=350, height=350) %>%
              plotly::layout(margin = list( l = 5, r = 5, b = 10, t = 10, pad = 0),
                             autosize = TRUE, 
                             plot_bgcolor  = "rgba(255, 255, 255, 1)", 
                             paper_bgcolor = "rgba(0, 0, 0, 0)") %>% 
              plotly::config(displayModeBar = FALSE) 
            
            sunplot_large <- plotly::plot_ly(sunplot, ids = ~ids, labels = ~labels, parents = ~parents, type = 'sunburst', width=900, height=900) %>%
              plotly::layout(margin = list( l = 20, r = 20, b = 20, t = 20, pad = 0),
                             autosize = TRUE, 
                             plot_bgcolor  = "rgba(255, 255, 255, 1)", 
                             paper_bgcolor = "rgba(0, 0, 0, 0)") %>% 
              plotly::config(displayModeBar = FALSE) 
            
            #-----------------------------------------------------------------------------------------------------------------Make other interactive plots
            if(length(unique(merged_input$rank_genus)) > 7){
              merged_input.g <- merged_input %>% group_by(rank_genus) %>% mutate(group_counts=n()) %>% ungroup() %>%
                arrange(desc(group_counts)) %>% mutate(rank_genus = ifelse(!rank_genus %in% unique(rank_genus)[1:7], "Others", rank_genus)) %>%
                arrange(desc(ID)) %>% mutate(group2 = ifelse(ID < species_threshold, "below", "above")) %>% mutate(bar_counts=rep(1)) %>%
                arrange(desc(group_counts))
              merged_input.g$rank_genus <- factor(merged_input.g$rank_genus, levels=c(unique(merged_input.g$rank_genus)))
              merged_input.g$filename <- factor(merged_input.g$filename, levels=c(unique(merged_input.g$filename)))
              #merged_input.g$ID <- factor(merged_input.g$ID, levels=c(unique(merged_input.g$ID)))
            } else {
              merged_input.g <- merged_input %>% group_by(rank_genus) %>% mutate(group_counts=n()) %>% ungroup() %>%
                arrange(desc(group_counts)) %>%
                arrange(desc(ID)) %>% mutate(group2 = ifelse(ID < species_threshold, "below", "above")) %>% mutate(bar_counts=rep(1)) %>%
                arrange(desc(group_counts))
              merged_input.g$rank_genus <- factor(merged_input.g$rank_genus, levels=c(unique(merged_input.g$rank_genus)))
              merged_input.g$filename <- factor(merged_input.g$filename, levels=c(unique(merged_input.g$filename)))
              #merged_input.g$ID <- factor(merged_input.g$ID, levels=c(unique(merged_input.g$ID)))
            }
            
            qual.range <- max(merged_input.g$phred_trim) - min(merged_input.g$phred_trim)
            if(qual.range == 0){phred_trim_adjust <- 0}
            if(qual.range != 0){phred_trim_adjust <- qual.range*0.01} #Jitter points by 1% of range
            
            g.scatter <- ggplot2::ggplot(merged_input.g, ggplot2::aes(x=ID , y=phred_trim, colour=ID)) + 
              ggplot2::geom_smooth(method='lm', formula= y~x, colour=scales::alpha("grey15", 0.3), linewidth=1.75, fill="#CCE6CC", alpha=0.7) + # 9F71A8 #40A499 #CCE6CC
              ggiraph::geom_point_interactive(ggplot2::aes(x=ID-0.005, y=phred_trim-phred_trim_adjust, group=ID, colour=ID, data_id = filename, tooltip=filename), 
                                              position=ggplot2::position_jitter(width=0, height=phred_trim_adjust),size = 0.7, hover_nearest = TRUE, alpha=0.5, colour="black") +
              #viridis::scale_colour_viridis(option="B", begin=0, end=0.8, direction=-1, alpha=0.5) +
              ggplot2::ylab("Phred quality score") + ggplot2::xlab("% identity to closest match") +
              ggplot2::theme(panel.background = ggplot2::element_blank(),
                             panel.border = ggplot2::element_rect(color ="black", fill = NA, linewidth = 0.5, linetype = 1),
                             panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                             legend.title = ggplot2::element_text(size=16),legend.text = ggplot2::element_text(size=12),
                             strip.background = ggplot2::element_blank(), strip.text.x = ggplot2::element_blank(),
                             axis.ticks.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_text(size=12), 
                             axis.title.y = ggplot2::element_text(size=16), axis.text.y=ggplot2::element_text(size=12),
                             axis.title.x= ggplot2::element_text(size=16)) + ggplot2::theme(legend.position="none")
            
            ggiraph::set_girafe_defaults(
              opts_hover_inv = ggiraph::opts_hover_inv(css = "opacity:0.2;"), 
              opts_hover = ggiraph::opts_hover(css = "stroke:orange;stroke-width:7;"),
              opts_tooltip = ggiraph::opts_tooltip(css = "padding:3px;background-color:#333333;color:white;"),
              opts_sizing = ggiraph::opts_sizing(rescale = FALSE),
              opts_toolbar = ggiraph::opts_toolbar(saveaspng = TRUE, pngname="isoQC", delay_mouseout = 5000, position="bottomleft", hidden=c("selection", "zoom"))
            )
            
            isoTAX.plots1 <- ggiraph::girafe(ggobj = g.scatter, width_svg = 3.5, height_svg = 4.75)
            
            g.id <- ggplot2::ggplot(merged_input.g, ggplot2::aes(x=rank_genus, y=bar_counts, fill=ID, colour=ID)) + 
              ggplot2::scale_x_discrete(limits=rev) +
              ggiraph::geom_bar_interactive(ggplot2::aes(fill=ID,
                                                         data_id=reorder(filename, ID),
                                                         tooltip = reorder(paste0(filename , "\n Closest type strain: ", closest_match, "\n % similarity: ", ID, "%"), ID)),
                                            stat="identity", width=0.9, colour=scales::alpha("grey85", 0), hover_nearest = FALSE) +
              #viridis::scale_fill_viridis(option="B", begin=0.2, direction=-1, alpha=0.6) +
              ggplot2::scale_fill_gradient(low="grey30", high="#CCE6CC") + # (option="B", begin=0.2, direction=-1, alpha=0.6) +
              ggplot2::scale_y_continuous(expand = c(0.02,0.02)) +
              ggplot2::coord_flip() + 
              ggplot2::ylab("Seqeunce counts per genus")+
              ggplot2::theme(panel.background = ggplot2::element_blank(),
                             panel.border = ggplot2::element_rect(color ="black", fill = NA, linewidth = 0.5, linetype = 1),
                             panel.grid.major = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(),
                             legend.title = ggplot2::element_text(size=16),legend.text = ggplot2::element_text(size=12),
                             strip.background = ggplot2::element_blank(), strip.text.x = ggplot2::element_blank(),
                             axis.ticks.x=ggplot2::element_blank(), axis.text.x=ggplot2::element_text(size=12), 
                             axis.title.y = ggplot2::element_blank(), axis.text.y=ggplot2::element_text(size=12, margin=ggplot2::margin(t = 0, r = 5, b = 0, l = 0, unit = "pt")),
                             axis.title.x= ggplot2::element_text(size=16))
            
            ggiraph::set_girafe_defaults(
              opts_hover_inv = ggiraph::opts_hover_inv(css = "opacity:0.2;"), 
              opts_hover = ggiraph::opts_hover(css = ''),
              opts_tooltip = ggiraph::opts_tooltip(css = "padding:3px;background-color:#333333;color:white;"),
              opts_sizing = ggiraph::opts_sizing(rescale = FALSE),
              opts_toolbar = ggiraph::opts_toolbar(saveaspng = TRUE, pngname="isoQC", delay_mouseout = 5000, position="topright", hidden=c("selection", "zoom"))
            )
            
            isoTAX.plots2 <- ggiraph::girafe(ggobj = g.id, width_svg = 9.5, height_svg = 4.75)
            
            #-----------------------------------------------------------------------------------------------------------------
            
            #Import S4 object
            uni.font <- 11
            
            html_output <- crosstalk::bscols(widths = 12,
                                             reactable(merged_input,
                                                       fullWidth=FALSE,
                                                       searchable = TRUE,
                                                       bordered = TRUE,
                                                       resizable =TRUE,
                                                       pageSizeOptions = c(10, 25, 50, 100, 1000),
                                                       defaultPageSize=15,
                                                       showPageSizeOptions = TRUE,
                                                       highlight = TRUE,
                                                       showSortable = TRUE,
                                                       compact=TRUE,
                                                       wrap = TRUE, #highlight=TRUE,
                                                       theme = reactableTheme(                  #borderColor = "#555",
                                                         headerStyle = list(
                                                           fontFamily = "sans-serif", fontWeight="bold", fontSize=uni.font,
                                                           "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                                                           "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)")),
                                                         style = list(
                                                           fontFamily = "sans-serif", fontWeight="normal",
                                                           "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                                                           "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"))),
                                                       style = list(minWidth = 1400),
                                                       #--------------------Universal settings
                                                       defaultColDef = colDef(vAlign = "center", style = list(whiteSpace = "nowrap", fontSize=uni.font, align="center", fontFamily = "sans-serif")),
                                                       #--------------------Column specific settings
                                                       columns = list(
                                                         warning = colDef(name="Warning", vAlign = "center", headerVAlign = "bottom", minWidth=170, #style = list(fontWeight="bold")),# cell = color_tiles(merged_input, color_ref = 'species_col')),
                                                                          style = function(value, index, alignment){
                                                                            if(!merged_input$warning[[index]]=="") {
                                                                              color.background <- grDevices::adjustcolor("#B20A2C", alpha.f=0.7)
                                                                              color.text <- "white"
                                                                            } else if(merged_input$warning[[index]]=="") {
                                                                              color.background <- NA
                                                                              color.text <- "black"  }
                                                                            list(background = color.background, color = color.text,
                                                                                 whiteSpace = "nowrap",
                                                                                 fontSize=uni.font-2,
                                                                                 align="center",
                                                                                 fontFamily = "sans-serif")}),
                                                         seqs_trim = colDef(name = "trimSeq",
                                                                            minWidth=70,
                                                                            style=list(whiteSpace = "nowrap",
                                                                                       fontSize=uni.font-2,
                                                                                       align="center",
                                                                                       fontFamily = "sans-serif")),
                                                         filename = colDef(name="Filename",
                                                                           minWidth= 150,
                                                                           style = function(value, index, alignment){
                                                                             if(!merged_input$warning[[index]]=="") {
                                                                               color.background <- grDevices::adjustcolor("#B20A2C", alpha.f=0.7)
                                                                               color.text <- "white"
                                                                             } else if(merged_input$warning[[index]]=="") {
                                                                               color.background <- NA
                                                                               color.text <- "black"  }
                                                                             list(background = color.background, color = color.text,
                                                                                  whiteSpace = "nowrap",
                                                                                  fontSize=uni.font-2,
                                                                                  align="center",
                                                                                  fontFamily = "sans-serif")}),
                                                         date = colDef(name = "Date",
                                                                       minWidth=75,
                                                                       style = function(value, index, alignment){
                                                                         if(!merged_input$warning[[index]]=="") {
                                                                           color.background <- grDevices::adjustcolor("#B20A2C", alpha.f=0.7)
                                                                           color.text <- "white"
                                                                         } else if(merged_input$warning[[index]]=="") {
                                                                           color.background <- NA
                                                                           color.text <- "black"  }
                                                                         list(background = color.background, color = color.text,
                                                                              whiteSpace = "nowrap",
                                                                              fontSize=uni.font-2,
                                                                              align="center",
                                                                              fontFamily = "sans-serif")}),
                                                         decision = colDef(name="Decision",  align="center", vAlign = "center",  minWidth=80,
                                                                           style = function(value, index, decision){
                                                                             if(merged_input$decision[[index]] == "Pass") {
                                                                               color.background <- grDevices::adjustcolor("#CCE6CC", alpha.f=1)
                                                                               color.text <- "black"
                                                                             } else if(!merged_input$decision[[index]] == "Pass") {
                                                                               color.background <- grDevices::adjustcolor("#B20A2C", alpha.f=0.7)
                                                                               color.text <- "white" }
                                                                             list(background = color.background, color = color.text,
                                                                                  whiteSpace = "nowrap",
                                                                                  fontSize=uni.font+4,
                                                                                  #fontWeight="bold",
                                                                                  align="center",
                                                                                  fontFamily = "sans-serif")}),
                                                         # Taxonomy columns
                                                         #-----------------------------------
                                                         rank_species = colDef(name="Species",
                                                                               vAlign = "center",
                                                                               headerVAlign = "top",
                                                                               minWidth=135,
                                                                               style = function(value, index, rank_species){
                                                                                 if(merged_input$ID[[index]] > species_threshold) {
                                                                                   color.background <- grDevices::adjustcolor("#CCE6CC", alpha.f=1)
                                                                                   color.text <- "black"
                                                                                 } else if(merged_input$ID[[index]] <= species_threshold){
                                                                                   color.background <- NA
                                                                                   color.text <- "black"
                                                                                 }
                                                                                 list(background = color.background, color = color.text,
                                                                                      whiteSpace = "nowrap",
                                                                                      fontSize=uni.font-2,
                                                                                      #fontWeight="bold",
                                                                                      align="center",
                                                                                      fontFamily = "sans-serif")}),
                                                         rank_genus = colDef(name="Genus",
                                                                             vAlign = "center",
                                                                             headerVAlign = "top",
                                                                             minWidth=75,
                                                                             style = function(value, index, rank_genus){
                                                                               if(merged_input$ID[[index]] > genus_threshold) {
                                                                                 color.background <- grDevices::adjustcolor("#CCE6CC", alpha.f=1)
                                                                                 color.text <- "black"
                                                                               } else if(merged_input$ID[[index]] <= genus_threshold){
                                                                                 color.background <- NA
                                                                                 color.text <- "black"
                                                                               }
                                                                               list(background = color.background, color = color.text,
                                                                                    whiteSpace = "nowrap",
                                                                                    fontSize=uni.font-2,
                                                                                    #fontWeight="bold",
                                                                                    align="center",
                                                                                    fontFamily = "sans-serif")}),
                                                         rank_family = colDef(name="Family",
                                                                              vAlign = "center",
                                                                              headerVAlign = "top",
                                                                              minWidth=75,
                                                                              style = function(value, index, rank_family){
                                                                                if(merged_input$ID[[index]] > family_threshold) {
                                                                                  color.background <- grDevices::adjustcolor("#CCE6CC", alpha.f=1)
                                                                                  color.text <- "black"
                                                                                } else if(merged_input$ID[[index]] <= family_threshold){
                                                                                  color.background <- NA
                                                                                  color.text <- "black"
                                                                                }
                                                                                list(background = color.background, color = color.text,
                                                                                     whiteSpace = "nowrap",
                                                                                     fontSize=uni.font-2,
                                                                                     #fontWeight="bold",
                                                                                     align="center",
                                                                                     fontFamily = "sans-serif")}),
                                                         rank_order = colDef(name="Order",
                                                                             vAlign = "center",
                                                                             headerVAlign = "top",
                                                                             minWidth=75,
                                                                             style = function(value, index, rank_order){
                                                                               if(merged_input$ID[[index]] > order_threshold) {
                                                                                 color.background <- grDevices::adjustcolor("#CCE6CC", alpha.f=1)
                                                                                 color.text <- "black"
                                                                               } else if(merged_input$ID[[index]] <= order_threshold){
                                                                                 color.background <- NA
                                                                                 color.text <- "black"
                                                                               }
                                                                               list(background = color.background, color = color.text,
                                                                                    whiteSpace = "nowrap",
                                                                                    fontSize=uni.font-2,
                                                                                    #fontWeight="bold",
                                                                                    align="center",
                                                                                    fontFamily = "sans-serif")}),
                                                         rank_class = colDef(name="Class",
                                                                             vAlign = "center",
                                                                             headerVAlign = "top",
                                                                             minWidth=75,
                                                                             style = function(value, index, rank_class){
                                                                               if(merged_input$ID[[index]] > class_threshold) {
                                                                                 color.background <- grDevices::adjustcolor("#CCE6CC", alpha.f=1)
                                                                                 color.text <- "black"
                                                                               } else if(merged_input$ID[[index]] <= class_threshold){
                                                                                 color.background <- NA
                                                                                 color.text <- "black"
                                                                               }
                                                                               list(background = color.background, color = color.text,
                                                                                    whiteSpace = "nowrap",
                                                                                    fontSize=uni.font-2,
                                                                                    #fontWeight="bold",
                                                                                    align="center",
                                                                                    fontFamily = "sans-serif")}),
                                                         rank_phylum = colDef(name="Phylum",
                                                                              vAlign = "center",
                                                                              headerVAlign = "top",
                                                                              minWidth=75,
                                                                              style = function(value, index, rank_phylum){
                                                                                if(merged_input$ID[[index]] > phylum_threshold) {
                                                                                  color.background <- grDevices::adjustcolor("#CCE6CC", alpha.f=1)
                                                                                  color.text <- "black"
                                                                                } else if(merged_input$ID[[index]] <= phylum_threshold){
                                                                                  color.background <- NA
                                                                                  color.text <- "black"
                                                                                }
                                                                                list(background = color.background, color = color.text,
                                                                                     whiteSpace = "nowrap",
                                                                                     fontSize=uni.font-2,
                                                                                     #fontWeight="bold",
                                                                                     align="center",
                                                                                     fontFamily = "sans-serif")}),
                                                         # Hidden columns
                                                         #-------------------
                                                         #warning = colDef(show=FALSE),
                                                         seq_warnings_col = colDef(show=FALSE),
                                                         species_colors = colDef(show=FALSE),
                                                         species_col = colDef(show=FALSE),
                                                         genus_col = colDef(show=FALSE),
                                                         family_col = colDef(show=FALSE),
                                                         class_col = colDef(show=FALSE),
                                                         order_col = colDef(show=FALSE),
                                                         phylum_col = colDef(show=FALSE),
                                                         Ns_raw = colDef(show=FALSE),
                                                         phred_raw = colDef(show=FALSE),
                                                         seqs_raw = colDef(show=FALSE),
                                                         length_raw = colDef(show=FALSE),
                                                         length_raw = colDef(show=FALSE),
                                                         phred_spark_raw = colDef(show=FALSE),
                                                         # Gauge widgets
                                                         #-------------------
                                                         closest_match = colDef(name = "Closest match (type strain)",
                                                                                minWidth=220,
                                                                                style=list(color = "black",
                                                                                           whiteSpace = "nowrap",
                                                                                           fontSize=uni.font-2,
                                                                                           align="center",
                                                                                           fontFamily = "sans-serif")),
                                                         ID = colDef(header= function(value, name) {
                                                           htmltools::div(title = "% sequence identity of query sequence to closest match sequence in database", value)},
                                                           align= "center",
                                                           vAlign = "center",
                                                           headerVAlign = "bottom",
                                                           width= 60,
                                                           cell = color_tiles(data=merged_input, colors="black", text_size=uni.font+4, bold_text=TRUE)),
                                                         NCBI_acc = colDef(width=80,
                                                                           style=list(color = "black",
                                                                                      whiteSpace = "nowrap",
                                                                                      fontSize=uni.font-2,
                                                                                      align="center",
                                                                                      fontFamily = "sans-serif")),
                                                         phred_trim = colDef(name = "trimQ",
                                                                             width= 60,
                                                                             headerVAlign = "top",
                                                                             cell = gauge_chart(merged_input,
                                                                                                fill_color = c('#e5f5e0', '#a1d99b', '#31a354'),
                                                                                                opacity = 0.5,
                                                                                                number_fmt = scales::comma,
                                                                                                text_size = 12,
                                                                                                max_value = max(c(merged_input$phred_raw,merged_input$phred_trim)),
                                                                                                show_min_max = FALSE)),
                                                         Ns_trim = colDef(name = "trimNs",
                                                                          width= 65,
                                                                          headerVAlign = "top",
                                                                          cell = gauge_chart(merged_input,
                                                                                             fill_color = c("white", "#B20A2C"), opacity = 0.5,
                                                                                             number_fmt = scales::comma,
                                                                                             #bold_text = TRUE,
                                                                                             text_size = 12,
                                                                                             max_value = max(c(merged_input$Ns_raw,merged_input$Ns_trim)),
                                                                                             show_min_max = FALSE)),
                                                         # Barplot widgets
                                                         #-------------------
                                                         length_trim = colDef(name = "trimLength",
                                                                              width=90,
                                                                              headerVAlign = "top",
                                                                              cell = data_bars(merged_input,
                                                                                               text_position = "inside-base", #inside-end / outside-base / above
                                                                                               text_size = 12,
                                                                                               #fill_color="#6f86ab",
                                                                                               fill_color = c('#FFF2D9','#FFE1A6','#FFCB66','#FFB627'),
                                                                                               fill_gradient = TRUE,
                                                                                               background = 'transparent',
                                                                                               number_fmt = scales::comma_format(accuracy = 1),
                                                                                               round_edges = FALSE, align_bars="left")) )) %>%  reactablefmtr::google_font(font_family = "Source Sans Pro") %>%
                                               htmlwidgets::prependContent(
                                                 shiny::tags$div(
                                                   shiny::tags$div(
                                                     shiny::tags$div(
                                                       "isoTAX output table", 
                                                       style = htmltools::css('font-size' = '24pt', 'font-weight' = 'bold',
                                                                              'text-align' = 'left', 'margin-bottom' = 0, 'padding-top' = '10px', 'padding-left' = '10px', 'vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("---------------------------------------------------", sep=""), 
                                                                     style = htmltools::css('font-size' = '14pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste('Project directory: "', unlist(strsplit(obj@input, '/'))[length(unlist(strsplit(obj@input, '/')))-1], '"', sep=""), 
                                                                     style = htmltools::css('font-size' = '10pt', 'font-weight' = 'bold', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("Date last updated: ", Sys.Date(), sep=""), 
                                                                     style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("---------------------------------------------------", sep=""), 
                                                                     style = htmltools::css('font-size' = '14pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("Settings:", sep=""), 
                                                                     style = htmltools::css('font-size' = '10pt', 'font-weight' = 'bold', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("quick_search = '", quick_search ,"'", sep=""),
                                                                     style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-top' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("db = '", db,"'", sep=""),
                                                                     style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("phylum_threshold = ", obj@phylum_threshold[1], sep=""),
                                                                     style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("class_threshold = ", obj@class_threshold[1], sep=""),
                                                                     style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("order_threshold = ",obj@order_threshold[1], sep=""),
                                                                     style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("family_threshold = ",obj@family_threshold[1], sep=""),
                                                                     style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("genus_threshold = ",obj@genus_threshold[1], sep=""),
                                                                     style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     shiny::tags$div(paste("species_threshold = ",obj@species_threshold[1], sep=""),
                                                                     style = htmltools::css('font-size' = '10pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                            'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '10px','vertical-align' = 'middle')),
                                                     style = htmltools::css(width = '28%')),
                                                   
                                                   shiny::tags$div(
                                                     shiny::tags$div(sunplot_small, style = htmltools::css(display='inline-flex', 'align' = 'left', 'margin-left' = '-22%',
                                                                                                           width='15%', 'vertical-align' = 'bottom'))
                                                   ),
                                                   shiny::tags$div(
                                                     shiny::tags$a(href = gsub("////", "///", paste("file:///",file.path(path, "lib", "sb_large.html"), sep="")), "Full Size", target="_blank"),
                                                     style = htmltools::css(display='inline-flex', width='10%', 'font-size' = '10pt', 'alpha' = '0.5', 'padding-top' = '15px','padding-left' = '10px',
                                                                            'font-weight' = 'normal','margin-left' = '-10%', 'z-index' = '1')
                                                   ),
                                                   shiny::tags$div(isoTAX.plots1, style = htmltools::css(display='inline-flex', width='15%', 'margin-left' = '-4%', 'margin-bottom' = '-10%')), 
                                                   shiny::tags$div(isoTAX.plots2, style = htmltools::css(display='inline-flex', width='20%', 'padding-left' = '50px', 'margin-bottom' = '-2%')),
                                                   style = htmltools::css(
                                                     width = '1400px',
                                                     display = 'inline-flex')))
            )

            sunplot_large1 <- shiny::tags$div(sunplot_large, 'align'= "center", 'height' = '100vh', 'margin'='0') 
            fname_html_sun <- file.path(path, "lib", "sb_large.html")
            htmltools::save_html(sunplot_large1, fname_html_sun)
            fname_html <- file.path(path, output)
            htmltools::save_html(html_output, fname_html)
            writeLines(gsub('<meta charset="utf-8"/>', paste('<meta charset="utf-8"/>\n<title>isoTAX > ', unlist(strsplit(obj@input, '/'))[length(unlist(strsplit(obj@input, '/')))-1],'</title>', sep=""), readLines(fname_html)), fname_html)
            pander::openFileInOS(fname_html)
            #return(html_output)
          })


#:::::::::::::::::::::::::::::::
# MAKE HTML #3: LIB TABLE
#:::::::::::::::::::::::::::::::
#' @export
#' @name export_html
#' @rdname export_html
#' @aliases export_html-isoLIB
#' @import crosstalk
#' @import htmltools
#' @importFrom crosstalk bscols
#' @importFrom crosstalk filter_checkbox
#' @importFrom crosstalk filter_select
#' @importFrom crosstalk filter_slider
#' @importFrom crosstalk SharedData
#' @importFrom pander openFileInOS
#' @importFrom dataui dui_sparkline
#' @importFrom dataui dui_sparkbandline
#' @importFrom dataui dui_sparklineseries
#' @importFrom htmltools save_html
#' @importFrom htmltools div
#' @importFrom htmltools browsable
#' @import reactable
#' @import reactablefmtr
#' @importFrom shiny tagList
#' @importFrom scales comma
#' @importFrom scales comma_format
#' @importFrom scales label_number


setMethod("export_html", "isoLIB",
          function(obj, method=NULL, group_cutoff=NULL){
            #Import checks
            if(is.null(method)){ message("The 'method' parameter is blank. Setting to 'NA' to enable graphing.")
              method <- 'NA'}
            if(is.null(group_cutoff)){ message("The 'group_cutoff' parameter is blank. Setting to 'NA' to enable graphing.")
              group_cutoff <- 'NA'}

            html_input2 <- S4_to_dataframe(obj) %>%
              mutate(phylum_col = ifelse(ID > phylum_threshold, "#31a354", "#FFFFFF00")) %>%
              mutate(class_col = ifelse(ID > class_threshold, "#31a354", "#FFFFFF00")) %>%
              mutate(order_col = ifelse(ID > order_threshold, "#31a354", "#FFFFFF00")) %>%
              mutate(family_col = ifelse(ID > family_threshold, "#31a354", "#FFFFFF00")) %>%
              mutate(genus_col = ifelse(ID > genus_threshold, "#31a354", "#FFFFFF00")) %>%
              mutate(species_col = ifelse(ID > species_threshold, "#31a354", "#FFFFFF00")) %>%
              dplyr::relocate(representative, .after=filename)
            
            #if(!grepl("maximum_likelihood|closest_species|dark_mode", method)){ method <- "NA"}
            #if(method=="closest_species"){ method <- "NA"}
            #if(is.null(group_cutoff)){ group_cutoff <- "NA"}
            
            html_input <- html_input2 %>% select(-input,
                                                 -phylum_threshold,
                                                 -class_threshold,
                                                 -order_threshold,
                                                 -family_threshold,
                                                 -genus_threshold,
                                                 -species_threshold,
                                                 -phylum_col,
                                                 -class_col,
                                                 -order_col,
                                                 -family_col,
                                                 -genus_col,
                                                 -species_col) %>% arrange(desc(representative))
            
            uni.font <- 9
            set.spacing <- "&nbsp;&nbsp;&nbsp;&nbsp; | &nbsp;&nbsp;&nbsp;&nbsp;"
            path <- getwd()
            output <- "03_isoLIB_results.html"
            
            data <- crosstalk::SharedData$new(html_input)
            
            
            if(length(unique(html_input2$ID)) <2  | length(unique(html_input2$length_trim)) <2){
              
              html_output <- crosstalk::bscols(widths = c(1,10),
                                               list(
                                                 #filter_checkbox("categories", "Categories", data, ~categories, inline = TRUE),
                                                 #crosstalk::filter_slider("length_trim", "Seq length", round=-1, ticks=FALSE, data, ~length_trim),
                                                 #crosstalk::filter_slider("ID", "% identity", data, round=2, ticks=FALSE, ~ID),
                                                 crosstalk::filter_checkbox("representative", "Sequence Group Representative (Rep)", data, ~representative, inline = FALSE),
                                                 crosstalk::filter_checkbox("date", "Date Sequenced", data, ~date, inline = FALSE)
                                               ),
                                               htmltools::browsable(
                                                 shiny::tagList(tags$button(
                                                   "Expand/collapse all",
                                                   onclick = "Reactable.toggleAllRowsExpanded('html_input')"
                                                 ),
                                                 shiny::tags$button("Download as CSV", onclick = "Reactable.downloadDataCSV('html_input')"),
                                                 #new_library_reactable <-
                                                 reactable(data,
                                                           fullWidth=FALSE,
                                                           searchable = TRUE,
                                                           bordered = TRUE,
                                                           resizable =TRUE,
                                                           pageSizeOptions = c(10, 20, 50, 100, 1000),
                                                           showPageSizeOptions = TRUE,
                                                           defaultPageSize=20,
                                                           highlight = TRUE,
                                                           showSortable = TRUE,
                                                           compact=TRUE,
                                                           style = list(minWidth = 1400),
                                                           theme = reactableTheme(
                                                             headerStyle = list(
                                                               fontFamily = "sans-serif",
                                                               fontWeight="bold",
                                                               fontSize=uni.font+1,
                                                               "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                                                               "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"))),
                                                           defaultColDef = colDef(minWidth= 75,
                                                                                  vAlign = "center",
                                                                                  style = list(whiteSpace = "nowrap",
                                                                                               fontSize=uni.font,
                                                                                               align="center",
                                                                                               fontFamily = "sans-serif")),
                                                           groupBy = "sequence_group",
                                                           elementId = "html_input",
                                                           #colDef options
                                                           #---------------
                                                           columns = list(
                                                             sequence_group = colDef(name="Sequence Group", minWidth=140,
                                                                                     grouped = JS("function(cellInfo) {
                                                                      if (cellInfo.subRows.length > 100) {
                                                                        return cellInfo.value + ' (' + cellInfo.subRows.length + ')'
                                                                      }
                                                                      return cellInfo.value
                                                                    }")),
                                                             #V10 = colDef(style = "white-space: nowrap;", aggregate = "unique"),
                                                             filename = colDef(name = "Filename",
                                                                               minWidth=120,
                                                                               aggregate = "count",
                                                                               style = list(whiteSpace = "nowrap",
                                                                                            fontSize=uni.font,
                                                                                            align="center",
                                                                                            fontFamily = "sans-serif")),
                                                             
                                                             representative = colDef(name = "Rep",
                                                                                     minWidth=70,
                                                                                     align="center"),
                                                             phred_trim = colDef(name = "Q",
                                                                                 minWidth= 35,
                                                                                 align="center",
                                                                                 aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                                          Ns_trim = colDef(name = "Ns",
                                                                           minWidth= 40,
                                                                           align="center",
                                                                           aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                                          length_trim = colDef(name = "Length",
                                                                               minWidth=70,
                                                                               align="left",
                                                                               aggregate = JS("function(values, rows){
                                                                                           return values[0]
                                                                      }"),
                                                                      cell = data_bars(html_input,
                                                                                       text_position = "inside-base",
                                                                                       text_size = 10,
                                                                                       fill_color = c('#FFF2D9','#FFE1A6','#FFCB66','#FFB627'),
                                                                                       fill_gradient = TRUE,
                                                                                       background = 'transparent',
                                                                                       number_fmt = scales::comma_format(accuracy = 1),
                                                                                       round_edges = FALSE, align_bars="left")),
                                                          seqs_trim = colDef(name="Sequence",
                                                                             minWidth= 75,
                                                                             style = list(whiteSpace = "nowrap",
                                                                                          fontSize=uni.font,
                                                                                          align="center",
                                                                                          fontFamily = "sans-serif"),
                                                                             aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                                          closest_match = colDef(name="Closest match (type strain)",
                                                                                 minWidth=220, #aggregate = "unique"
                                                                                 aggregate = JS("function(values, rows){
                                                             return values[0]
                                                            }")),
                                                          date = colDef(name= "Date",
                                                                        minWidth= 58,
                                                                        style = list(whiteSpace = "nowrap",
                                                                                     fontSize=uni.font,
                                                                                     align="center",
                                                                                     fontFamily = "sans-serif"),
                                                                        aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                                          NCBI_acc = colDef(minWidth=80,
                                                                            aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                                          ID = colDef(align= "center",
                                                                      vAlign = "center",
                                                                      headerVAlign = "top",
                                                                      minWidth= 50,
                                                                      aggregate = JS("function(values, rows){
                                                             return values[0]
                                                            }"),
                                                                      cell = color_tiles(data=html_input,
                                                                                         colors="black",
                                                                                         number_fmt = scales::label_number(accuracy = 0.1),
                                                                                         text_size=uni.font+3,
                                                                                         bold_text=TRUE),
                                                                      header = function(value, name) {
                                                                        htmltools::div(title = "% sequence identity of query sequence to closest match sequence in database", value)}),
                                                          rank_phylum = colDef(name="Phylum",
                                                                               align = "left",
                                                                               vAlign = "center",
                                                                               headerVAlign = "top",
                                                                               minWidth=110,
                                                                               style = list(whiteSpace = "nowrap",
                                                                                            fontSize=uni.font,
                                                                                            align="left",
                                                                                            fontFamily = "sans-serif"),
                                                                               aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="phylum_col")),
                                                          rank_class = colDef(name="Class",
                                                                              align = "left",
                                                                              vAlign = "center",
                                                                              headerVAlign = "top",
                                                                              minWidth=110,
                                                                              style = list(whiteSpace = "nowrap",
                                                                                           fontSize=uni.font,
                                                                                           align="left",
                                                                                           fontFamily = "sans-serif"),
                                                                              aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="class_col")),
                                                          rank_order = colDef(name="Order",
                                                                              align = "left",
                                                                              vAlign = "center",
                                                                              headerVAlign = "top",
                                                                              minWidth=110,
                                                                              aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="order_col"),
                                                              style = list(whiteSpace = "nowrap",
                                                                           fontSize=uni.font,
                                                                           align="left",
                                                                           fontFamily = "sans-serif")
                                                          ),
                                                          rank_family = colDef(name="Family",
                                                                               align = "left",
                                                                               vAlign = "center",
                                                                               headerVAlign = "top",
                                                                               minWidth=110,style = list(whiteSpace = "nowrap",
                                                                                                         fontSize=uni.font,
                                                                                                         align="left",
                                                                                                         fontFamily = "sans-serif"),
                                                                               aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="family_col")),
                                                          rank_genus = colDef(name="Genus",
                                                                              align = "left",
                                                                              vAlign = "center",
                                                                              headerVAlign = "top",
                                                                              minWidth=110, style = list(whiteSpace = "nowrap",
                                                                                                         fontSize=uni.font,
                                                                                                         align="left",
                                                                                                         fontFamily = "sans-serif"),
                                                                              aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="genus_col")),
                                                          rank_species = colDef(name="Species",
                                                                                align = "left",
                                                                                vAlign = "center",
                                                                                headerVAlign = "top",
                                                                                minWidth=175,
                                                                                style = list(whiteSpace = "nowrap",
                                                                                             fontSize=uni.font,
                                                                                             align="left",
                                                                                             fontFamily = "sans-serif"),
                                                                                aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="species_col"))
                                                           )) %>% reactablefmtr::google_font(font_family = "Source Sans Pro") %>%
                                                   htmlwidgets::prependContent(
                                                     shiny::tags$div(
                                                       shiny::tags$div(
                                                         shiny::tags$div(paste('Project directory: "', unlist(strsplit(obj@input, '/'))[length(unlist(strsplit(obj@input, '/')))-1], '"', sep=""), 
                                                                         style = htmltools::css('font-size' = '12pt', 'font-weight' = 'bold', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                                'margin-bottom' = 0, 'padding-top' = '10px', 'padding-left' = '0px','vertical-align' = 'middle')),
                                                         shiny::tags$div(paste("Date last updated: ", Sys.Date(), sep=""), 
                                                                         style = htmltools::css('font-size' = '12pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                                'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '0px','vertical-align' = 'middle')),
                                                         shiny::tags$div(paste("method = '", method, "' | group_cutoff = '", group_cutoff, "'", sep=""),
                                                                         style = htmltools::css('font-size' = '12pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                                'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '0px','vertical-align' = 'middle'))),
                                                       shiny::tags$div("isoLIB output table", 
                                                                       style = htmltools::css('font-size' = '24pt', 'font-weight' = 'bold',
                                                                                              'text-align' = 'left', 'margin-bottom' = 0, 'padding-top' = '10px', 'padding-left' = '28%', 'vertical-align' = 'middle')),
                                                       style = htmltools::css(
                                                         width = '1400px',
                                                         display = 'flex'))) %>%
                                                   add_subtitle("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------", font_size = 14, font_style="normal", font_weight="normal") %>%
                                                   add_subtitle(paste("Total sequences: " ,length(obj@filename) , set.spacing,
                                                                      "Total strain groups: ", length(obj@representative[obj@representative == "yes"]) , set.spacing,
                                                                      "Mean Phred quality: ", format(round(mean(obj@phred_trim), 0), 2),  set.spacing,
                                                                      "Mean no. Ns: ", format(round(mean(obj@Ns_trim), 2), 2),  set.spacing,
                                                                      "Mean length: ", round(mean(obj@length_trim), 0),  set.spacing,
                                                                      "No. unique phyla: ", round(length(unique(obj@rank_phylum)), 0),  set.spacing,
                                                                      "No. unique classes: ", round(length(unique(obj@rank_class)), 0),  set.spacing,
                                                                      "No. unique orders: ", round(length(unique(obj@rank_order)), 0),  set.spacing,
                                                                      "No. unique families: ", round(length(unique(obj@rank_family)), 0),  set.spacing,
                                                                      "No. unique genera: ", round(length(unique(obj@rank_genus)), 0),  set.spacing,
                                                                      "No. unique species: ", round(length(unique(obj@rank_species)), 0),
                                                                      sep=""), font_size = 14, font_style="normal", font_weight="normal", align="left") %>%
                                                   add_subtitle("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------", font_size = 14, font_style="normal", font_weight="normal")
                                                 ,
                                                 )))
            } else {html_output <- crosstalk::bscols(widths = c(1,10),
                                                     list(
                                                       #filter_checkbox("categories", "Categories", data, ~categories, inline = TRUE),
                                                       crosstalk::filter_slider("length_trim", "Sequence length", round=-1, ticks=FALSE, data, ~length_trim),
                                                       crosstalk::filter_slider("ID", "% identity", data, round=2, ticks=FALSE, ~ID),
                                                       crosstalk::filter_checkbox("representative", "Sequence Group Representative (Rep)", data, ~representative, inline = FALSE),
                                                       crosstalk::filter_checkbox("date", "Date Sequenced", data, ~date, inline = FALSE)
                                                     ),
                                                     
                                                     htmltools::browsable(
                                                       shiny::tagList(
                                                         shiny::tags$button(
                                                           "Expand/collapse all",
                                                           onclick = "Reactable.toggleAllRowsExpanded('html_input')"
                                                         ),
                                                         shiny::tags$button("Download as CSV", onclick = "Reactable.downloadDataCSV('html_input')"),
                                                         #new_library_reactable <-
                                                         reactable(data,
                                                                   fullWidth=FALSE,
                                                                   searchable = TRUE,
                                                                   bordered = TRUE,
                                                                   resizable =TRUE,
                                                                   pageSizeOptions = c(10, 20, 50, 100, 1000),
                                                                   showPageSizeOptions = TRUE,
                                                                   defaultPageSize=20,
                                                                   highlight = TRUE,
                                                                   showSortable = TRUE,
                                                                   compact=TRUE,
                                                                   style = list(minWidth = 1400),
                                                                   theme = reactableTheme(
                                                                     headerStyle = list(
                                                                       fontFamily = "sans-serif",
                                                                       fontWeight="bold",
                                                                       fontSize=uni.font+1,
                                                                       "&:hover[aria-sort]" = list(background = "hsl(0, 0%, 96%)"),
                                                                       "&[aria-sort='ascending'], &[aria-sort='descending']" = list(background = "hsl(0, 0%, 96%)"))),
                                                                   defaultColDef = colDef(minWidth= 75,
                                                                                          vAlign = "center",
                                                                                          style = list(whiteSpace = "nowrap",
                                                                                                       fontSize=uni.font,
                                                                                                       align="center",
                                                                                                       fontFamily = "sans-serif")),
                                                                   groupBy = "sequence_group",
                                                                   elementId = "html_input",
                                                                   #colDef options
                                                                   #---------------
                                                                   columns = list(
                                                                     sequence_group = colDef(name="Sequence Group", minWidth=140,
                                                                                             grouped = JS("function(cellInfo) {
                                                                      if (cellInfo.subRows.length > 100) {
                                                                        return cellInfo.value + ' (' + cellInfo.subRows.length + ')'
                                                                      }
                                                                      return cellInfo.value
                                                                    }")),
                                                                    #V10 = colDef(style = "white-space: nowrap;", aggregate = "unique"),
                                                                    filename = colDef(name = "Filename",
                                                                                      minWidth=120,
                                                                                      aggregate = "count",
                                                                                      style = list(whiteSpace = "nowrap",
                                                                                                   fontSize=uni.font,
                                                                                                   align="center",
                                                                                                   fontFamily = "sans-serif")),
                                                                    representative = colDef(name = "Rep",
                                                                                            minWidth=70,
                                                                                            align="center"),
                                                                    phred_trim = colDef(name = "Q",
                                                                                        minWidth= 35,
                                                                                        align="center",
                                                                                        aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                                          Ns_trim = colDef(name = "Ns",
                                                                           minWidth= 40,
                                                                           align="center",
                                                                           aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                                          length_trim = colDef(name = "Length",
                                                                               minWidth=70,
                                                                               align="left",
                                                                               aggregate = JS("function(values, rows){
                                                                                           return values[0]
                                                                      }"),
                                                                      cell = data_bars(html_input,
                                                                                       text_position = "inside-base",
                                                                                       text_size = 10,
                                                                                       fill_color = c('#FFF2D9','#FFE1A6','#FFCB66','#FFB627'),
                                                                                       fill_gradient = TRUE,
                                                                                       background = 'transparent',
                                                                                       number_fmt = scales::comma_format(accuracy = 1),
                                                                                       round_edges = FALSE, align_bars="left")),
                                                          seqs_trim = colDef(name="Sequence",
                                                                             minWidth= 75,
                                                                             style = list(whiteSpace = "nowrap",
                                                                                          fontSize=uni.font,
                                                                                          align="center",
                                                                                          fontFamily = "sans-serif"),
                                                                             aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                                          closest_match = colDef(name="Closest match (type strain)",
                                                                                 minWidth=220, #aggregate = "unique"
                                                                                 aggregate = JS("function(values, rows){
                                                             return values[0]
                                                            }")),
                                                          date = colDef(name= "Date",
                                                                        minWidth= 58,
                                                                        style = list(whiteSpace = "nowrap",
                                                                                     fontSize=uni.font,
                                                                                     align="center",
                                                                                     fontFamily = "sans-serif"),
                                                                        aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                                          NCBI_acc = colDef(minWidth=80,
                                                                            aggregate = JS("function(values, rows){
                                                                               return values[0]
                                                          }")),
                                                          ID = colDef(align= "center",
                                                                      vAlign = "center",
                                                                      headerVAlign = "top",
                                                                      minWidth= 50,
                                                                      aggregate = JS("function(values, rows){
                                                             return values[0]
                                                            }"),
                                                                      cell = color_tiles(data=html_input,
                                                                                         colors="black",
                                                                                         number_fmt = scales::label_number(accuracy = 0.1),
                                                                                         text_size=uni.font+3,
                                                                                         bold_text=TRUE),
                                                                      header = function(value, name) {
                                                                        htmltools::div(title = "% sequence identity of query sequence to closest match sequence in database", value)}),
                                                          rank_phylum = colDef(name="Phylum",
                                                                               align = "left",
                                                                               vAlign = "center",
                                                                               headerVAlign = "top",
                                                                               minWidth=110,
                                                                               style = list(whiteSpace = "nowrap",
                                                                                            fontSize=uni.font,
                                                                                            align="left",
                                                                                            fontFamily = "sans-serif"),
                                                                               aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="phylum_col")),
                                                          rank_class = colDef(name="Class",
                                                                              align = "left",
                                                                              vAlign = "center",
                                                                              headerVAlign = "top",
                                                                              minWidth=110,
                                                                              style = list(whiteSpace = "nowrap",
                                                                                           fontSize=uni.font,
                                                                                           align="left",
                                                                                           fontFamily = "sans-serif"),
                                                                              aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="class_col")),
                                                          rank_order = colDef(name="Order",
                                                                              align = "left",
                                                                              vAlign = "center",
                                                                              headerVAlign = "top",
                                                                              minWidth=110,
                                                                              aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="order_col"),
                                                              style = list(whiteSpace = "nowrap",
                                                                           fontSize=uni.font,
                                                                           align="left",
                                                                           fontFamily = "sans-serif")
                                                          ),
                                                          rank_family = colDef(name="Family",
                                                                               align = "left",
                                                                               vAlign = "center",
                                                                               headerVAlign = "top",
                                                                               minWidth=110,style = list(whiteSpace = "nowrap",
                                                                                                         fontSize=uni.font,
                                                                                                         align="left",
                                                                                                         fontFamily = "sans-serif"),
                                                                               aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="family_col")),
                                                          rank_genus = colDef(name="Genus",
                                                                              align = "left",
                                                                              vAlign = "center",
                                                                              headerVAlign = "top",
                                                                              minWidth=110, style = list(whiteSpace = "nowrap",
                                                                                                         fontSize=uni.font,
                                                                                                         align="left",
                                                                                                         fontFamily = "sans-serif"),
                                                                              aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="genus_col")),
                                                          rank_species = colDef(name="Species",
                                                                                align = "left",
                                                                                vAlign = "center",
                                                                                headerVAlign = "top",
                                                                                minWidth=175,
                                                                                style = list(whiteSpace = "nowrap",
                                                                                             fontSize=uni.font,
                                                                                             align="left",
                                                                                             fontFamily = "sans-serif"),
                                                                                aggregate = JS("function(values, rows){
                                                              return values[0]
                                                              }"),
                                                              cell = color_tiles(html_input2, color_ref="species_col"))
                                                                   )) %>% reactablefmtr::google_font(font_family = "Source Sans Pro") %>%
                                                           
                                                           htmlwidgets::prependContent(
                                                               shiny::tags$div(
                                                                 shiny::tags$div(
                                                                   shiny::tags$div(paste('Project directory: "', unlist(strsplit(obj@input, '/'))[length(unlist(strsplit(obj@input, '/')))-1], '"', sep=""), 
                                                                                   style = htmltools::css('font-size' = '12pt', 'font-weight' = 'bold', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                                          'margin-bottom' = 0, 'padding-top' = '10px', 'padding-left' = '0px','vertical-align' = 'middle')),
                                                                   shiny::tags$div(paste("Date last updated: ", Sys.Date(), sep=""), 
                                                                                   style = htmltools::css('font-size' = '12pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                                          'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '0px','vertical-align' = 'middle')),
                                                                   shiny::tags$div(paste("method = '", method, "' | group_cutoff = '", group_cutoff, "'", sep=""),
                                                                                   style = htmltools::css('font-size' = '12pt', 'font-weight' = 'normal', 'font-style' = 'normal', 'text-align' = 'left',
                                                                                                          'margin-bottom' = 0, 'padding-top' = '0px', 'padding-left' = '0px','vertical-align' = 'middle'))),
                                                                 shiny::tags$div("isoLIB output table", 
                                                                   style = htmltools::css('font-size' = '24pt', 'font-weight' = 'bold',
                                                                                          'text-align' = 'left', 'margin-bottom' = 0, 'padding-top' = '10px', 'padding-left' = '28%', 'vertical-align' = 'middle')),
                                                                 style = htmltools::css(
                                                                   width = '1400px',
                                                                   display = 'flex'))) %>%
                                                           add_subtitle("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------", font_size = 14, font_style="normal", font_weight="normal") %>%
                                                           add_subtitle(paste("Total sequences: " ,length(obj@filename) , set.spacing,
                                                                              "Total strain groups: ", length(obj@representative[obj@representative == "yes"]) , set.spacing,
                                                                              "Mean Phred quality: ", format(round(mean(obj@phred_trim), 0), 2),  set.spacing,
                                                                              "Mean no. Ns: ", format(round(mean(obj@Ns_trim), 2), 2),  set.spacing,
                                                                              "Mean length: ", round(mean(obj@length_trim), 0),  set.spacing,
                                                                              "No. unique phyla: ", round(length(unique(obj@rank_phylum)), 0),  set.spacing,
                                                                              "No. unique classes: ", round(length(unique(obj@rank_class)), 0),  set.spacing,
                                                                              "No. unique orders: ", round(length(unique(obj@rank_order)), 0),  set.spacing,
                                                                              "No. unique families: ", round(length(unique(obj@rank_family)), 0),  set.spacing,
                                                                              "No. unique genera: ", round(length(unique(obj@rank_genus)), 0),  set.spacing,
                                                                              "No. unique species: ", round(length(unique(obj@rank_species)), 0),
                                                                              sep=""), font_size = 14, font_style="normal", font_weight="normal", align="left") %>%
                                                           add_subtitle("---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------", font_size = 14, font_style="normal", font_weight="normal")
                                                         ,
                                                       )))
            }
            fname_html <- file.path(path, output)
            htmltools::save_html(html_output, fname_html)
            writeLines(gsub('<meta charset="utf-8"/>', paste('<meta charset="utf-8"/>\n<title>isoLIB > ', unlist(strsplit(obj@input, '/'))[length(unlist(strsplit(obj@input, '/')))-1],'</title>', sep=""), readLines(fname_html)), fname_html)
            #writeLines(gsub('<meta charset="utf-8"/>', '<meta charset="utf-8"/>\n<title>isoLIB</title>', readLines(fname_html)), fname_html)
            writeLines(gsub('amp;', '', readLines(fname_html)), fname_html)
            writeLines(gsub('">---', ';white-space:nowrap">---', readLines(fname_html)), fname_html)
            writeLines(gsub('">Total seq', ';white-space:nowrap">Total seq', readLines(fname_html)), fname_html)
            pander::openFileInOS(fname_html)
            #return(html_output)
          })






