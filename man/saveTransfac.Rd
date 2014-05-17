\name{saveTransfac}
\alias{saveTransfac}
\title{Save a PWM matrix object to transfac format}
\usage{
saveTransfac(pwm, file.out = tempfile("transfac"), type = "aa")
}
\arguments{
  \item{pwm}{PWM matrix object}

  \item{file.out}{where the transfac matrix is written}

  \item{type}{'aa', 'dna' or 'rna' depending on the
  namespace}
}
\description{
This saves an already generated PWM matrix object in R to
transfac format, which can be read in by RWebLogo
}
\keyword{transfac}

