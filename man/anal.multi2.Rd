\name{anal.multi2}
\alias{anal.multi2}

\title{Simulations to test QTL analyses}

\description{Run multiple simulations of a backcross, analyzed using
  ANOVA, CIM, forward selection (with both BIC and permutation
  tests) and MCMC.}

\usage{
anal.multi2(n.sim=1000, cim.steps=c(3,5,7,9,11), bic.mult=c(2,2.5,3),
            max.steps=13, n.perm=1000, alpha=0.05, thresh=NULL, drop=1.5,
            n.ind=100, n.mar=rep(11,9), mar.sp=rep(10,10*9),
            n.mcmc=1000, mcmc.bic=2.56,
            qtl.chr=c(1,1,2,2,3,4,5), qtl.mar=c(4,8,4,8,6,4,1),
            qtl.dist=rep(0,7), qtl.eff=c(1,1,1,-1,1,1,1)/2,
            sigma=1)
}
\arguments{
  \item{n.sim}{Number of simulation replicates.}
  \item{cim.steps}{Number of steps in forward selection prior to CIM analysis.}
  \item{bic.mult}{Multipliers for BIC.}
  \item{max.steps}{Maximum number of steps in forward selection.}
  \item{n.perm}{Number of permutation replicates in permutation tests.}
  \item{alpha}{Significance level in permutation tests.}
  \item{thresh}{LOD thresholds for ANOVA and CIM.  Should be length 1 or
    1+length(cim.steps).  If NULL and the chromosome structure and
    number of individuals match simulations that I've already done, we
    use my estimated LOD thresholds.}
  \item{drop}{Drop in LOD that is required between inferred QTLs.
    Should be length 1 or 1+length(cim.steps).}
  \item{n.ind}{Number of backcross individuals.}
  \item{n.mar}{Vector indicating the number of markers on each chromosome.}
  \item{mar.sp}{Vector of length sum(n.mar)-n.chr-1, giving the
    inter-marker spacings in cM.}
  \item{n.mcmc}{Number of steps in MCMC}
  \item{mcmc.bic}{Delta value in BIC in the MCMC runs}
  \item{qtl.chr}{Chromosomes on which qtls are sitting.}
  \item{qtl.mar}{Markers to the left of each qtl, numbered 1, 2, ...}
  \item{qtl.dist}{Distance between qtl and marker to left, in cM.}
  \item{qtl.eff}{QTL effects.}
  \item{sigma}{Residual (environmental) standard deviation.}
}

\value{
  A list of length 3+n.cim+n.bic, whose components correspond to the
  results of ANOVA, CIM, forward selection with BIC, forward
  selection with permutation tests, and MCMC with BIC.  Each component
  of this list is itself a list of length n.sim, with each component
  giving the chromosome numbers and marker numbers of the identified
  QTLs.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\seealso{ \code{\link{simbc}}, \code{\link{sim.null}},
  \code{\link{sim.mcmc}}, \code{\link{anal.multi}} }

\examples{
\dontrun{results <- anal.multi2(n.sim=10)}
}

\references{
  Broman, K. W. (1997) Identifying quantitative trait loci in experimental
  crosses.  PhD dissertation, Department of Statistics, University of
  California, Berkeley.
}

\keyword{models}
\keyword{regression}
\keyword{datagen}
