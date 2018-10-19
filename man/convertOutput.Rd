\name{convertOutput}
\alias{convertOutput}

\title{Convert output of anal.leaps, anal.mcmc or perm}

\description{Convert the output of anal.leaps, anal.mcmc or perm for use
  with the function which.correct.}

\usage{
convertOutput(object, n.mar=rep(11,9))
}

\arguments{
  \item{object}{An object of the form output from
    \code{\link{anal.leaps}}, \code{\link{anal.mcmc}},
    \code{\link{perm}}; namely, a vector of marker identifiers,
    as integers, (or a list whose first component is such a vector).}
  \item{n.mar}{Vector indicating the number of markers on each chromosome.}
}

\value{
  A matrix whose columns are identified QTLs, and rows are chromosome
  and marker numbers, respectively.
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\seealso{ \code{\link{anal.leaps}},
  \code{\link{anal.mcmc}}, \code{\link{perm}},
  \code{\link{which.correct}} }

\examples{
dat <- simbc(n.ind=250)
res1 <- anal.leaps(dat, "forward", 2.1, 27)
res2 <- anal.leaps(dat, "backward", 2.1, 27)
res <- lapply(list(res1,res2), convertOutput)
out <- which.correct(res)
}

\keyword{manip}
