% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoALL.R
\name{isoALL}
\alias{isoALL}
\title{Perform all commands in one step.}
\usage{
isoALL(
  input = NULL,
  export_html = TRUE,
  export_csv = TRUE,
  export_fasta = TRUE,
  export_fasta_revcomp = FALSE,
  quick_search = FALSE,
  db = "16S",
  iddef = 2,
  phylum_cutoff = 75,
  class_cutoff = 78.5,
  order_cutoff = 82,
  family_cutoff = 86.5,
  genus_cutoff = 96.5,
  species_cutoff = 98.7,
  include_warnings = FALSE,
  strain_group_cutoff = 0.995,
  merge = FALSE
)
}
\arguments{
\item{input}{Directory path(s) containing .ab1 files. If more than one, provivde as list (e.g. 'input=c("/path/to/directory1","/path/to/directory2")')}

\item{export_html}{(Default=TRUE) Output the results as an HTML file}

\item{export_csv}{(Default=TRUE) Output the results as a CSV file.}

\item{export_fasta}{(Default=TRUE) Output the sequences in a FASTA file.}

\item{export_fasta_revcomp}{(Default=FALSE) Output the sequences in reverse complement form in a fasta file. This is useful in cases where sequencing was done using the reverse primer and thus the orientation of input sequences needs reversing.}

\item{quick_search}{(Default=FALSE) Whether or not to perform a comprehensive database search (i.e. optimal global alignment).
If TRUE, performs quick search equivalent to setting VSEARCH parameters "--maxaccepts 100 --maxrejects 100".
If FALSE, performs comprehensive search equivalent to setting VSEARCH parameters "--maxaccepts 0 --maxrejects 0"}

\item{db}{(Default="16S") Select database option(s) including "16S" (for searching against the NCBI Refseq targeted loci 16S rRNA database),
"ITS" (for searching against the NCBI Refseq targeted loci ITS  database. For combined databases in cases where input sequences are dervied from
bacteria and fungi, select "16S|ITS".}

\item{iddef}{Set pairwise identity definition as per VSEARCH definitions (Default=2, and is recommended for highest taxonomic accuracy)
(0) CD-HIT definition: (matching columns) / (shortest sequence length).
(1) Edit distance: (matching columns) / (alignment length).
(2) Edit distance excluding terminal gaps (default definition).
(3) Marine Biological Lab definition counting each gap opening (internal or terminal) as a single mismatch, whether or not the gap was extended: 1.0- ((mismatches + gap openings)/(longest sequence length)).
(4) BLAST definition, equivalent to --iddef 1 for global pairwise alignments.}

\item{phylum_cutoff}{Percent cutoff for phylum rank demarcation}

\item{class_cutoff}{Percent cutoff for class rank demarcation}

\item{order_cutoff}{Percent cutoff for order rank demarcation}

\item{family_cutoff}{Percent cutoff for family rank demarcation}

\item{genus_cutoff}{Percent cutoff for genus rank demarcation}

\item{species_cutoff}{Percent cutoff for species rank demarcation}

\item{include_warnings}{(Default=FALSE) Whether or not to keep sequences with poor alignment warnings from Step 2 'isoTAX' function. Set TRUE to keep warning sequences, and FALSE to remove warning sequences.}

\item{strain_group_cutoff}{(Default=0.995) Similarity cutoff (0-1) for delineating between strain groups. 1 = 100\% identical/0.995=0.5\% difference/0.95=5.0\% difference/etc.}

\item{merge}{If TRUE, combines isoLIB output files consecutively in the order they are listed. Default=FALSE performs all the steps (isoQC/isoTAX/isoLIB) on each directory separately.}

\item{verbose}{(Default=FALSE) Output progress while script is running.}

\item{files_manual}{(Default=NULL) For testing purposes only. Specify a list of files to run  as filenames without extensions, rather than the whole directory format. Primarily used for testing, use at your own risk.}

\item{exclude}{(Default=NULL) For testing purposes only. Excludes files of interest from input directory.}

\item{min_phred_score}{(Default=20) Do not accept trimmed sequences with a mean Phred score below this cutoff}

\item{min_length}{(Default=200) Do not accept trimmed sequences with sequence length below this number}

\item{sliding_window_cutoff}{(Default=NULL) Quality trimming parameter (M2) for wrapping SangerRead function in sangeranalyseR package. If NULL, implements auto cutoff for Phred score (recommended), otherwise set between 1-60.}

\item{sliding_window_size}{(Default=15) Quality trimming parameter (M2) for wrapping SangerRead function in sangeranalyseR package. Recommended range between 5-30.}

\item{date}{Set date "YYYY_MM_DD" format. If NULL, attempts to parse date from .ab1 file}
}
\value{
Returns a list of \code{\link{class-isoLIB}} class objects.
}
\description{
This function effectively wraps isoQC, isoTAX, and isoLIB steps into a single command for convenience. Input can be a single directory
or a list of directories to process at once. If multiple directories are provided, the resultant libraries can be sequentially merged together
by toggling the parameter 'merge=TRUE'. All other respective parameters from the wrapped functions can be passed through this command.
. The The respective input parameters from
the wrappred  can be passed through this command with exception of the .creates a strain library by grouping closely related strains of interest based on sequence similarity.
For adding new sequences to an already-established strain library, specify the .CSV file path of the older strain library using the 'old_lib_csv" parameter.
}
\examples{
#Set path to directory containing example .ab1 files
fpath1 <- system.file("extdata/abif_examples/rocket_salad", package = "isolateR")

#Run isoALL function with default settings
isoALL(input=fpath1)
}
\seealso{
\code{\link{isoQC}}, \code{\link{isoTAX}}, \code{\link{isoLIB}}
}