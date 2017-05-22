\name{plot_res}
\alias{plot_res}

\title{Plot simulation results}

\description{Plot a summary of the results of
  \code{\link[qtlsim]{anal.multi}} or \code{\link[qtlsim]{sim.mcmc}}.}

\usage{
plot_res(x, \dots, ylim=NULL)
}

\arguments{
  \item{x}{Output of \code{\link[qtlsim]{which.correct}}.}
  \item{\dots}{Arguments to be passed to \code{\link[graphics]{par}}.}
  \item{ylim}{Y-axis limits.  If NULL, (0,1) is used.}
}

\value{
  A plot of the proportion of times each QTL was identified by each
  method, and the proportion of linked and unlinked extraneous loci.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{anal.multi}}, \code{\link{anal.mcmc}},
  \code{\link{which.correct}} }

\examples{
\dontrun{output <- anal.multi()}
\dontrun{tab <- which.correct(output)}
\dontrun{plot_res(tab)}
}

\keyword{hplot}
