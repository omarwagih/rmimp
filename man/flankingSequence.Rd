\name{flankingSequence}
\alias{flankingSequence}
\title{Get flanking sequences of a position.}
\usage{
flankingSequence(seqs, inds, flank = 7, empty.char = "-")
}
\arguments{
  \item{seqs}{Charachter vector of sequences. If only one
  sequence is provided, indices from \code{inds} are
  assumed to all be from the same sequence.}

  \item{inds}{Numerical vector of positions corresponding
  to the sequences provided in \code{seqs}.}

  \item{flank}{Value indicating the number of charachters
  to extract, before and after an index}

  \item{empty.char}{Charachter used to replace out of bound
  flanking sequences}
}
\description{
This function obtains the flanking sequence at one or more
position. Out of bound indices are replaced by a blank
charachter.
}
\examples{
# One sequence and one index. Central charachter is 'B'
flankingSequence(seqs='ABC', inds=2, flank=1)
# An example showing the use of empty.char
flankingSequence(seqs='ABC', inds=2, flank=5)
# An example with muliple sequences and indicies
flankingSequence(seqs=c('ABC', 'XYZ'), inds=c(2, 1), flank=1)
}
\keyword{flank}
\keyword{sequence}

