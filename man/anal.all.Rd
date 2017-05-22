\name{anal.all}
\alias{anal.all}

\title{Analyze a backcross}

\description{Analyze a backcross using ANOVA, CIM, and forward selection.}

\usage{
anal.all(dat, cim.steps=7, max.steps=20)
}

\arguments{
  \item{dat}{The data: a list with components \code{geno} (the genotype
    data, as a matrix) and \code{pheno} (the phenotype data, as a vector).}
  \item{cim.steps}{Number of steps in forward selection prior to CIM analysis.}
  \item{max.steps}{Maximum number of steps in forward selection.}
}

\value{
  A list with three components: LOD scores from ANOVA, LOD scores from
  CIM, and a matrix of marker indices and RSS from forward selection
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{simbc}}, \code{\link{anal.leaps}},
  \code{\link{anal.mcmc}}, \code{\link{perm}} }


\examples{
dat <- simbc()
results <- anal.all(dat)
}

\references{
  Broman, K. W. (1997) Identifying quantitative trait loci in experimental
  crosses.  PhD dissertation, Department of Statistics, University of
  California, Berkeley.
}

\keyword{models}
\keyword{regression}
