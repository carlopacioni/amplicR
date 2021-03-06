% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.proc.R
\name{data.proc}
\alias{data.proc}
\title{Data processing}
\usage{
data.proc(
  dir.in = NULL,
  dir.out = NULL,
  bp = 0,
  truncQ = 2,
  qrep = FALSE,
  dada = TRUE,
  pool = FALSE,
  plot.err = FALSE,
  chim = TRUE,
  orderBy = "abundance",
  verbose = TRUE
)
}
\arguments{
\item{dir.in}{The directory where the fastq files are located. If NULL 
(default) an interactive window is used to select a folder}

\item{dir.out}{The directory where to save the results. If NULL (default) 
then \code{dir.out <- dir.in}}

\item{bp}{An integer indicating the expected length (base-pairs) of the
reads. If zero (default) no truncation is applied}

\item{truncQ}{Truncate reads at the first instance of a quality score less 
than or equal to truncQ when conducting quality filtering. See
\code{\link[dada2]{fastqFilter}} for details}

\item{qrep}{Logical. Should the quality report be generated? (default 
\code{FALSE})}

\item{dada}{Logical. Should the dada analysis be conducted? (default 
\code{TRUE})}

\item{pool}{Logical. Should samples be pooled together prior to sample 
inference? (default \code{FALSE}). See \code{\link[dada2]{dada}} for
details}

\item{plot.err}{Logical. Whether error rates obtained from \code{dada} should
be plotted}

\item{chim}{Logical. Should the bimera search and removal be performed? 
(default \code{TRUE})}

\item{orderBy}{Character vector specifying how the returned sequence table 
should be sorted. Default "abundance". See 
\code{\link[dada2]{makeSequenceTable}} for details}

\item{verbose}{Logical. Whether information on progress should be outputted 
(default: TRUE)}
}
\value{
Return a list with several elements:
  
  \itemize{ \item $luniseqsFinal: A list with unique sequences (names) that
  were retained at completion of \code{data.proc} and their abundance
  (values). \item $lsummary: A list where each element is a summary of the
  number of reads that were retained in each step. This can be converted in a
  \code{data.frame} by running the following
  
  \code{summary <- plyr::join_all(lsummary, by="Sample", type="left")} \item
  $stable: The sequence table \item $seq_list: The sequences and matching 
  sequence IDs \item $call: The function call}
  
  Several files are also returned, these include: \itemize{ \item
  Seq_table.csv Sequence table \item Seq_list.csv List of sequences and their
  matching IDs \item data.proc.summary.csv A summary of the number of reads
  that were retained in each step \item data.proc.rda R data file containing
  the list returned by data.proc (see above) }
}
\description{
\code{data.proc} is a function to process (quality checking, error and 
chimeras filtering) data from a NGS run after these have been deconvoluted.
}
\details{
\code{data.proc} locates the .fastq files (can be compressed) in the 
directory indicated in \code{dir.in}. If the directory path is not provided, 
this will be selected using an interactive window.

It is currently limited to single-reads and assumes that adapters, primers 
and indexes have been already removed and that each file represents a sample.

The \code{data.proc} pipeline is as follows: fastq files are read in. A 
filter is applied to truncate reads at the first instance of a quality score 
less than \code{truncQ}, remove reads  that are of low quality (currently the
threshold is hard-coded and reads are discarded if the expected errors is
higher than 3 - from documentation in the R package \code{dada2}, the
expected errors are calculated from the nominal definition of the quality
score: EE = sum(10^(-Q/10)) - and remove reads that (after truncation) do not
match the target length. A quality report can be (optionally) generated with the
R package \code{ShortReads} to verify the quality of the reads retained
after this step. Reads are then dereplicated. Optionally, the dada (Callahan
et al 2015) algoritm is applied and bimeras are searched and removed with
default settings of the relative functions in the package \code{dada2}. The
sequences that were retained at completion of \code{data.proc} are saved in
fasta files in the subfolder "Final_seqs" and a .csv with a summary of the
number of reads that have been retained in each step is also written. These
two outputs are also returned at the end of the function.

The sequence data handling is done by using functionalities from the packages
\code{dada2} and \code{ShortRead}, so make sure to cite them (in addition to 
\code{amplicR} of course!) if you report your results in a paper or report.
}
\examples{
# Select the directory where the example data are stored
example.data <- system.file("extdata", "HTJ", package="amplicR")
# Select a temporary directory where to store the outputs
out <- tempdir()

HTJ.test <- data.proc(example.data, out, bp=140)

# To clean up the temp directory
unlink(out, recursive=TRUE)
}
\references{
Benjamin J Callahan, Paul J McMurdie, Michael J Rosen, Andrew W 
  Han, Amy J Johnson, Susan P Holmes (2015). DADA2: High resolution sample 
  inference from amplicon data.
}
\seealso{
\code{\link[dada2]{dada}}, \code{\link[dada2]{makeSequenceTable}}
}
