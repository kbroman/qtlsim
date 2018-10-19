\name{combine.sims}
\alias{combine.sims}

\title{Combine multiple simulation results}

\description{Combine output of multiple runs of
  \code{\link{anal.multi}}.}

\usage{
combine.sims(\dots)
}

\arguments{
  \item{\dots}{Two or more objects, of the form output from
  \code{\link{anal.multi}}.}
}

\value{
  A single object, of the form output from
  \code{\link{anal.multi}}.
}

\author{Karl W Broman, \email{broman@wisc.edu}}

\seealso{ \code{\link{anal.multi}},
  \code{\link{combine.mcmc}} }

\examples{
\dontrun{output1 <- anal.multi(n.sim=10)}
\dontrun{output2 <- anal.multi(n.sim=10)}
\dontrun{output <- combine.sims(output1, output2)}
}

\keyword{manip}
