\name{sim.null}
\alias{sim.null}

\title{Simulate under the null hypothesis of no QTLs.}

\description{Simulates genotype and phenotype data for a backcross
  experiment with no QTLs, and applies ANOVA and CIM, in order to
  estimate LOD thresholds.}

\usage{
sim.null(n.ind=100, n.mar=rep(11,9), mar.sp=rep(10,10*9),
         cim.steps=c(3,5,7,9,11), n.sim=10000)
}

\arguments{
  \item{n.ind}{Number of backcross individuals.}
  \item{n.mar}{Vector indicating the number of markers on each chromosome.}
  \item{mar.sp}{Vector of length sum(n.mar)-n.chr-1, giving the
    inter-marker spacings in cM.}
  \item{cim.steps}{Vector giving number of steps of forward selection to
    perform prior to CIM.}
  \item{n.sim}{Number of simulation replicates to perform.}
}

\value{
  A matrix of size [n.sim x (1+length(cim.steps))] giving the maximum
  LOD score for ANOVA and CIM (with the different values of cim.steps)
  for each simulation replicate.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{simbc}}, \code{\link{anal.multi}} }

\examples{
\dontrun{output <- sim.null()}
}

\references{
  Broman, K. W. (1997) Identifying quantitative trait loci in experimental
  crosses.  PhD dissertation, Department of Statistics, University of
  California, Berkeley.
}

\keyword{datagen}
\keyword{models}
\keyword{regression}
