\name{sim.mcmc}
\alias{sim.mcmc}

\title{Simulations to test QTL analyses}

\description{Run multiple simulations of a backcross, analyzed using
  several model selection techniques, including MCMC.}

\usage{
sim.mcmc(n.sim=1000, max.steps=28, bic.mult=2.5, n.steps=1000,
         n.ind=100, n.mar=rep(11,9), mar.sp=rep(10,10*9),
         qtl.chr=c(1,1,2,2,3,4,5), qtl.mar=c(4,8,4,8,6,4,1),
         qtl.dist=rep(0,7), qtl.eff=c(1,1,1,-1,1,1,1)/2,
         sigma=1, only.two=FALSE)
}

\arguments{
  \item{n.sim}{Number of simulation replicates.}
  \item{max.steps}{Maximum number of steps in forward/backward selection.}
  \item{bic.mult}{Multiplier for BIC.}
  \item{n.steps}{Number of steps in MCMC.}
  \item{n.ind}{Number of backcross individuals.}
  \item{n.mar}{Vector indicating the number of markers on each chromosome.}
  \item{mar.sp}{Vector of length sum(n.mar)-n.chr-1, giving the
    inter-marker spacings in cM.}
  \item{qtl.chr}{Chromosomes on which qtls are sitting.}
  \item{qtl.mar}{Markers to the left of each qtl, numbered 1, 2, ...}
  \item{qtl.dist}{Distance between qtl and marker to left, in cM.}
  \item{qtl.eff}{QTL effects.}
  \item{sigma}{Residual (environmental) standard deviation.}
  \item{only.two}{If TRUE, only do forward selection and MCMC}
}

\value{
  A list of three components.  The first is a list of length 5, whose
  components correspond to the 
  results of forward and backward selection, the same followed by
  branch-and-bound, with BIC, and MCMC.  (Each component of this list is
  itself a list of length n.sim, with each component giving the
  chromosome numbers and marker numbers of the identified QTLs.)
  The second component of the main list gives the BIC values
  corresponding to the identified models.  The last component gives the
  step at which the final model from MCMC was first observed, plus the
  number of distinct models observed.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} \cr
\url{http://www.biostat.wisc.edu/~kbroman/software} }

\seealso{ \code{\link{simbc}}, \code{\link{anal.mcmc}}, 
         \code{\link{anal.leaps}},
	 \code{\link{anal.multi}} }

\examples{
\dontrun{results <- sim.mcmc(n.sim=10)}
}
                    
\references{
  Broman, K. W. (1997) Identifying quantitative trait loci in experimental
  crosses.  PhD dissertation, Department of Statistics, University of
  California, Berkeley.
}

\keyword{models}
\keyword{regression}
\keyword{datagen}
