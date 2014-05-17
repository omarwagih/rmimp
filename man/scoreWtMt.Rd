\name{scoreWtMt}
\alias{scoreWtMt}
\title{Score wt and mt sequences for a pwm}
\usage{
scoreWtMt(pwm, mut_ps, is.kinase.pwm = T, thresh.bg = 1, thresh.fg = 0,
  thresh.log2 = 0)
}
\arguments{
  \item{pwm}{Position weight matrix of interest}

  \item{mut_ps}{psnvs data frame containing wt and mt
  sequences computed from pSNVs function}

  \item{is.kinase.pwm}{TRUE if pwm is that of a kinase}

  \item{thresh.bg}{Anything below this threshold is
  considered a negative hit}

  \item{thresh.fg}{Anything above this threshold is
  considered a positive hit}

  \item{thresh.log2}{Threshold for the absolute value of
  log ratio. Anything less than this value is discarded.}
}
\description{
Score wt and mt sequences for a pwm
}
\keyword{mut}
\keyword{psites}
\keyword{score}
\keyword{snp}

