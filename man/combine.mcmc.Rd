\name{combine.mcmc}
\alias{combine.mcmc}

\title{Combine multiple simulation results}

\description{Combine output of multiple runs of
  \code{\link{sim.mcmc}}.}

\usage{
combine.mcmc(\dots)
}

\arguments{
  \item{\dots}{Two or more objects, of the form output from
  \code{\link{sim.mcmc}}.}
}

\value{
  A single object, of the form output from
  \code{\link{sim.mcmc}}.
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\seealso{ \code{\link{sim.mcmc}},
  \code{\link{combine.sims}} }

\examples{
\dontrun{output1 <- sim.mcmc(n.sim=10)}
\dontrun{output2 <- sim.mcmc(n.sim=10)}
\dontrun{output <- combine.mcmc(output1, output2)}
}

\keyword{manip}
