\name{perm}
\alias{perm}

\title{Forward selection w/ permutation test}

\description{Analyze a backcross using forward selection with
  permutation tests to choose the size of the model.}

\usage{
perm(dat, n.perm=1000, alpha=0.05)
}

\arguments{
  \item{dat}{The data: a list with components \code{geno} (the genotype
    data, as a matrix) and \code{pheno} (the phenotype data, as a vector).}
  \item{n.perm}{Number of permutation replicates to perform at each state.}
  \item{alpha}{Level of significance to use in each permutation test.}
}

\value{
  A vector of marker indices.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} \cr
\url{http://www.biostat.wisc.edu/~kbroman/software} }

\seealso{ \code{\link{anal.all}}, \code{\link{anal.leaps}},
         \code{\link{anal.mcmc}}, \code{\link{simbc}} }

\examples{
\dontrun{dat <- simbc()}
\dontrun{results <- perm(dat)}
}
                    
\references{
  Broman, K. W. (1997) Identifying quantitative trait loci in experimental
  crosses.  PhD dissertation, Department of Statistics, University of
  California, Berkeley.
}


\keyword{models}
\keyword{regression}
