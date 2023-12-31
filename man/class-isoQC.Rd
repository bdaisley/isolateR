% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_class_methods.R
\docType{class}
\name{class-isoQC}
\alias{class-isoQC}
\title{isoQC Class Object}
\value{
Returns an class-isoQC object.
}
\description{
S4 wrapper for \code{\link{isoQC}} function. Access data via S4 slot functions.
}
\section{Slots}{

\describe{
\item{\code{date}}{Character string containing run date from each of the input Sanger sequence .ab1 files ("YYYY_MM_DD" format).}

\item{\code{filename}}{Character string containing input filenames.}

\item{\code{trim.start.pos}}{Numeric string containing trimming position start point.}

\item{\code{trim.end.pos}}{Numeric string containing trimming position end point.}

\item{\code{phred_spark_raw}}{List containing per nucleotide Phred score values for each sequence}

\item{\code{phred_raw}}{Numeric string containing mean Phred scores before trimming.}

\item{\code{phred_trim}}{Numeric string containing mean Phred scores after trimming.}

\item{\code{Ns_raw}}{Numeric string containing count of N's before trimming.}

\item{\code{Ns_trim}}{Numeric string containing count of N's after trimming.}

\item{\code{length_raw}}{Numeric string containing sequence length before trimming.}

\item{\code{length_trim}}{Numeric string containing sequence length after trimming.}

\item{\code{seqs_raw}}{Character string containing sequences before trimming.}

\item{\code{seqs_trim}}{Character string containing sequence after trimming.}

\item{\code{decision}}{Character string containing decision (PASS/FAIL) information based on \code{\link{isoQC}} 'min_phred_score' and 'min_length cutoffs'.}

\item{\code{input}}{Character string containing input directory information.}
}}

\seealso{
\code{\link{isoQC}}
}
