\name{PWM}
\alias{PWM}
\title{Construct position weight matrix}
\usage{
PWM(seqs, pseudocount = 0.001, relative.freq = T, type = "AA",
  priors = AA_PRIORS_HUMAN)
}
\arguments{
  \item{seqs}{Aligned sequences all of the same length}

  \item{pseudocount}{Pseudocount factor. Final pseudocount
  is background probability * this factor}

  \item{relative.freq}{TRUE if each column should be
  divided by the sum}

  \item{type}{Type of sequences 'AA' or 'DNA'}

  \item{priors}{Named character vector containing priors of
  amino acids.}
}
\description{
Makes a position weight matrix given aligned sequences.
}
\examples{
# No examples
}
\keyword{construct}
\keyword{pwm}

