\name{mss}
\alias{mss}
\title{Compute matrix similarity score as described in MATCH algorithm}
\usage{
mss(seqs, pwm, is.kinase.pwm = T, na.rm = F, ignore.ind = 8)
}
\arguments{
  \item{seqs}{Sequences to be scored}

  \item{pwm}{Position weight matrix}

  \item{is.kinase.pwm}{TRUE if PWM is that of a kinase}

  \item{na.rm}{Remove NA scores?}
}
\description{
Computes matrix similarity score of a PWM with a k-mer.
Score ranges from 0-1, as described in [PMID: 12824369]
}
\examples{
# No Examples
}
\keyword{match}
\keyword{mss}
\keyword{pwm}
\keyword{tfbs}

