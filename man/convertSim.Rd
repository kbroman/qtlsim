\name{convertSim}
\alias{convertSim}

\title{Convert simulated data into form for R/qtl}

\description{Convert data simulated with simbc() into
  the form used by R/qtl.}

\usage{
convertSim(object, n.chr=9, chr.len=100)
}

\arguments{
  \item{object}{A list with two components: the genotype data as a matrix,
    and the phenotype data as a vector.}
  \item{n.chr}{Number of chromosomes in the data.}
  \item{chr.len}{Length of the chromosomes, in cM.}
}

\value{
  An object of class \code{cross}.  See \code{\link[qtl]{read.cross}} for
  details.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu}}

\seealso{\code{\link[qtl]{read.cross}}, \code{\link{simbc}} }

\examples{
dat <- simbc()
cross <- convertSim(dat, 9, chr.len=100)
}

\keyword{manip}
