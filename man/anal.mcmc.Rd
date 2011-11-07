\name{anal.mcmc}
\alias{anal.mcmc}

\title{Analyze a backcross by MCMC}

\description{Analyze a backcross using an MCMC model selection method.}

\usage{
anal.mcmc(dat, bic.mult=2.5, n.steps=1000, start=NULL)
}

\arguments{
  \item{dat}{The data: a list with components \code{geno} (the genotype
    data, as a matrix) and \code{pheno} (the phenotype data, as a vector).}
  \item{bic.mult}{Multiplier for BIC, for choosing the size of the model.}
  \item{n.steps}{Number of steps to take in the Markov chain.}
  \item{start}{Model at which to start the chain: a vector giving the
    columns of the genotype matrix to use as an initial state.}
}

\value{
  A list giving the marker numbers chosen, the BIC criterion for the
  chosen model, the iteration of MCMC at which this model was first
  seen, and vectors giving the number of QTLs in the model and the BIC
  criterion, at each step of the chain
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} \cr
\url{http://www.biostat.wisc.edu/~kbroman/software} }

\seealso{ \code{\link{simbc}}, \code{\link{anal.all}},
         \code{\link{anal.leaps}}, \code{\link{perm}} }

\examples{
dat <- simbc()
results <- anal.mcmc(dat)
}
                    
\references{
  Broman, K. W. (1997) Identifying quantitative trait loci in experimental
  crosses.  PhD dissertation, Department of Statistics, University of
  California, Berkeley.
}

\keyword{models}
\keyword{regression}
