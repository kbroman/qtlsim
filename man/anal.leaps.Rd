\name{anal.leaps}
\alias{anal.leaps}

\title{Analyze a backcross by model selection methods}

\description{Analyze a backcross using forward or backward selection
             or either of these followed by branch-and-bound.}

\usage{
anal.leaps(dat, method=c("forward","backward","forwback",
                         "forwardleap","backwardleap"),
           bic.mult=2.5, max.steps=30)
}

\arguments{
  \item{dat}{The data: a list with components \code{geno} (the genotype
    data, as a matrix) and \code{pheno} (the phenotype data, as a vector).}
  \item{method}{Model search approach to use: forward selection,
                backward elimination, forward selection followed by
        backward elimination, or forward or backward selection
        followed by branch-and-bound.}
  \item{bic.mult}{Multiplier for BIC, for choosing the size of the model.}
  \item{max.steps}{Maximum number of steps in forward selection.}
}

\value{
  A list whose first component is a vector indicating the marker columns
  selected and second component is the corresponding BIC-delta value.
}

\author{Karl W Broman, \email{broman@wisc.edu}}


\seealso{ \code{\link{simbc}}, \code{\link{anal.all}},
         \code{\link{anal.mcmc}}, \code{\link{perm}} }

\examples{
dat <- simbc(n.ind=100,qtl.eff=c(1,1,1,-1,1,1,1)*0.75)
results <- anal.leaps(dat,method="forwback",bic.mult=2.56)
}

\references{
  Broman, K. W. (1997) Identifying quantitative trait loci in experimental
  crosses.  PhD dissertation, Department of Statistics, University of
  California, Berkeley.
}

\keyword{models}
\keyword{regression}
