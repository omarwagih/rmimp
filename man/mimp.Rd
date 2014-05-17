\name{mimp}
\alias{mimp}
\title{Predict the impact of single variants on phosphorylation.}
\usage{
mimp(muts, seqs, psites, perc.bg = 90, perc.fg = 10, thresh.log2 = 0,
  display.results = T)
}
\arguments{
  \item{muts}{Mutation data file: a space delimited text
  file or data frame containing two columns (1) gene and
  (1) mutation. Example: \tabular{ll}{ TP53 \tab R282W\cr
  CTNNB1 \tab S33C\cr CTNNB1 \tab S37F\cr }}

  \item{seqs}{Sequence data file containing protein
  sequences in FASTA format OR named list of sequences
  where each list element is the uppercase sequence and the
  name of each element is that of the protein. Example:
  list(TP53="ABCXYZ", CDK2="HJKEWR")}

  \item{psites}{Phosphorylation data file (optional): a
  space delimited text file containing positions of
  phosphorylation sites. Example: \tabular{ll}{ TP53 \tab
  280\cr CTNNB1 \tab 29\cr CTNNB1 \tab 44\cr }}

  \item{perc.bg}{Percentile value between 0 - 100. This
  value is used to compute a threshold, β from the negative
  (background) distribution of scores.  By default this is
  the 90th percentile of the background distribution of
  scores. Anything below the threshold is considered a
  negative hit.}

  \item{perc.fg}{Percentile value between 0 - 100. This
  value is used to compute a threshold, α from the positive
  (foreground) distribution of scores.  By default this is
  the 10th percentile of the foreground distribution of
  scores. Anything above the threshold is considered a
  positive hit.}

  \item{thresh.log2}{Threshold for the absolute value of
  log ratio. Anything less than this value is discarded
  (default: 0).}
}
\value{
The data is returned in a \code{data.frame} with the
following columns: \item{gene}{gene with the rewiring
event} \item{mut}{mutation causing the rewiring event}
\item{psite_pos}{position of the central residue of the
phosphosite} \item{mut_dist}{distance of the mutation from
the central phosphosite} \item{wt}{sequence of the wildtype
phosphosite (before the mutation) } \item{mt}{sequence of
the mutated phosphosite (after the mutation)}
\item{score_wt}{matrix similarity score of the wildtype
phosphosite } \item{score_mt}{matrix similarity score of
the mutated phosphosite } \item{log_ratio}{Log2 ratio
between mutant and wildtype scores. A high positive log
ratio represents a high confidence gain-of-signaling event.
A high negative log ratio represents a high confidence
loss-of-signaling event.} \item{pwm}{name of the kinase
being rewiried} \item{perc_wt}{Percentile rank of the wt
score} \item{perc_mt}{Percentile rank of the mutant score}
}
\description{
This function takes in mutation, sequence and
phosphorylation data to predict the impact the mutation has
on phosphorylation.
}
\examples{
# Get the path to example mutation data
mut.file = system.file("extdata", "mutation_data.txt", package = "MIMP")

# Get the path to example FASTA sequence data
seq.file = system.file("extdata", "sequence_data.txt", package = "MIMP")

# View the files in a text editor
browseURL(mut.file)
browseURL(seq.file)

# Run rewiring analysis
results = mimp(mut.file, seq.file, display.results=TRUE)

# Show head of results
head(results)
}
\keyword{kinase}
\keyword{mimp}
\keyword{mutation}
\keyword{phosphorylation}
\keyword{psites}
\keyword{pwm}
\keyword{rewiring}
\keyword{snp}
\keyword{snv}

