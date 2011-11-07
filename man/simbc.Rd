\name{simbc}
\alias{simbc}

\title{Simulate a backcross}

\description{Simulates genotype and phenotype data for a backcross
  experiment with QTLs exactly at marker postions, and acting
  additively.}

\usage{
simbc(n.ind=100, n.mar=rep(11,9), mar.sp=rep(10,10*9),
      qtl.chr=c(1,1,2,2,3,4,5), qtl.mar=c(4,8,4,8,6,4,1),
      qtl.dist=rep(0,7), qtl.eff=c(1,1,1,-1,1,1,1)/2,
      sigma=1)
}

\arguments{
  \item{n.ind}{Number of backcross individuals.}
  \item{n.mar}{Vector indicating the number of markers on each chromosome.}
  \item{mar.sp}{Vector of length sum(n.mar)-n.chr-1, giving the
    inter-marker spacings in cM.}
  \item{qtl.chr}{Chromosomes on which qtls are sitting.}
  \item{qtl.mar}{Markers to the left of each qtl, numbered 1, 2, ...}
  \item{qtl.dist}{Distance between qtl and marker to left, in cM.}
  \item{qtl.eff}{QTL effects.}
  \item{sigma}{Residual (environmental) standard deviation.}
}

\value{
  A list with two components: \code{geno}, a matrix containing the
  genotype data, and \code{pheno}, a vector of phenotypes.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} \cr
\url{http://www.biostat.wisc.edu/~kbroman/software} }

\seealso{ \code{\link{convertSim}} }

\examples{
dat <- simbc()
}
                    
\keyword{datagen}
