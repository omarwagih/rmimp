\name{pSNVs}
\alias{pSNVs}
\title{Find phosphorylation related variants (pSNVs)}
\usage{
pSNVs(muts, psites, seqs, flank = 7, multicore = F)
}
\arguments{
  \item{muts}{Mutation data as data frame of two columns
  (1) name of gene or protein (2) mutation in the format
  X123Y, where X is the reference amino acid and Y is the
  alternative amino acid.}

  \item{psites}{Phosphorylation data as a data frame of two
  columns (1) name of gene or protein (2) Position of the
  phosphorylated residue}

  \item{seqs}{Sequence data as a name list. Names of the
  list correspond to the gene or protein name. Each entry
  contains the collapsed sequence.}

  \item{flank}{Number of amino acids flanking the psite to
  be considered}

  \item{multicore}{If true, will use mclapply to speed
  things up!}
}
\description{
Given mutation data and psites, find variants that exist in
the flanking regions of the psite
}
\examples{
# No examples
}
\keyword{mutation}
\keyword{psites}
\keyword{psnv}
\keyword{snp}

