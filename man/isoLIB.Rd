% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isoLIB_function.R
\name{isoLIB}
\alias{isoLIB}
\title{Generate new strain library or add to existing one.}
\usage{
isoLIB(
  input = NULL,
  old_lib_csv = NULL,
  export_html = TRUE,
  export_csv = TRUE,
  include_warnings = FALSE,
  strain_group_cutoff = 0.995,
  phylum_cutoff = 75,
  class_cutoff = 78.5,
  order_cutoff = 82,
  family_cutoff = 86.5,
  genus_cutoff = 96.5,
  species_cutoff = 98.7
)
}
\arguments{
\item{input}{Path of CSV output file from isoTAX step.}

\item{old_lib_csv}{Optional: Path of CSV output isoLIB file or combined isoLIB file from previous run(s)}

\item{export_html}{(Default=TRUE) Output the results as an HTML file}

\item{export_csv}{(Default=TRUE) Output the results as a CSV file.}

\item{include_warnings}{(Default=FALSE) Whether or not to keep sequences with poor alignment warnings from Step 2 'isoTAX' function. Set TRUE to keep warning sequences, and FALSE to remove warning sequences.}

\item{strain_group_cutoff}{(Default=0.995) Similarity cutoff (0-1) for delineating between strain groups. 1 = 100\% identical/0.995=0.5\% difference/0.95=5.0\% difference/etc.}

\item{phylum_cutoff}{Percent cutoff for phylum rank demarcation}

\item{class_cutoff}{Percent cutoff for class rank demarcation}

\item{order_cutoff}{Percent cutoff for order rank demarcation}

\item{family_cutoff}{Percent cutoff for family rank demarcation}

\item{genus_cutoff}{Percent cutoff for genus rank demarcation}

\item{species_cutoff}{Percent cutoff for species rank demarcation}
}
\value{
Returns an isoLIB class object. Default taxonomic cutoffs for phylum (75.0), class (78.5), order (82.0), family (86.5), genus (96.5), and species (98.7) demarcation are based on Yarza et al. 2014, Nature Reviews Microbiology (DOI:10.1038/nrmicro3330)
}
\description{
This function creates a strain library by grouping closely related strains of interest based on sequence similarity.
For adding new sequences to an already-established strain library, specify the .CSV file path of the older strain library using the 'old_lib_csv" parameter.
}
\examples{
#Set path to directory containing example .ab1 files
fpath1 <- system.file("extdata/abif_examples/rocket_salad", package = "isolateR")

#Step 1: Run isoQC function with default settings
isoQC.S4 <- isoQC(input=fpath1)

#Step 2: Run isoTAX function with default settings
fpath2 <- file.path(fpath1, "isolateR_output/01_isoQC_trimmed_sequences_PASS.csv")
isoTAX.S4 <- isoTAX(input=fpath2)

#Step 3: Run isoLIB function with default settings
fpath3 <- file.path(fpath1, "isolateR_output/02_isoTAX_results.csv")
isoLIB.S4 <- isoLIB(input=fpath3)

#Show summary statistics
isoLIB.S4
}
\seealso{
\code{\link{isoTAX}}, \code{\link{isoLIB}}
}
