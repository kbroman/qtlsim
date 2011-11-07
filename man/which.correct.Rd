\name{which.correct}
\alias{which.correct}

\title{Determine which QTLs were correctly identified in an analysis}

\description{Given the results of a model selection procedure and the
  true underlying locations of QTLs, determine which QTLs were correctly
  identified and how many were extraneous.}

\usage{
which.correct(results, truth=rbind(c(1,1,2,2,3,4,5),c(4,8,4,8,6,4,1)),
              within=1)
}

\arguments{
  \item{results}{A list (or a list of lists), with each component in the
    inner lists being either \code{numeric(0)} or a matrix whose columns
    are identified QTLs, and rows are chromosome and marker numbers, 
    respectively.} 
  \item{truth}{Matrix whose columns are the true QTLs, with rows giving
    the chromosome and marker numbers, respectively.}
  \item{within}{Distance (in number of markers) allowable between true
    and inferred QTLs for the inferred QTL to be considered correct.}
}

\value{
  A list of matrices, with the matrices of size n.sim x (n.qtl+2).
  The first set of columns in each matrix indicate which QTLs were
  correctly identified.  The last two columns give the number of
  extraneous loci on the same chromosome as a QTL and the number
  of unlinked extraneous loci.
}

\author{Karl W Broman, \email{kbroman@biostat.wisc.edu} \cr
\url{http://www.biostat.wisc.edu/~kbroman/software} }

\seealso{ \code{\link{anal.multi}},
  \code{\link{anal.mcmc}} }

\examples{
\dontrun{dat <- simbc()}
\dontrun{output1 <- anal.multi(dat)}
\dontrun{output2 <- which.correct(output1)}
}
                    
\keyword{manip}

