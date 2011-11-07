######################################################################
#
# sims.R: Simulations of a backcross
#
# copyright (c) 2001-3, Karl W Broman
# July, 2001; June, 2003
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/qtlsim package
#
# simbc:       simulate backcross data 
#
######################################################################

simbc <-
function(n.ind=100, n.mar = rep(11,9), mar.sp = rep(10,10*9),
         qtl.chr = c(1,1,2,2,3,4,5), qtl.mar = c(4,8,4,8,6,4,1),
         qtl.dist = rep(0,7), qtl.eff = c(1,1,1,-1,1,1,1)/2,
         sigma=1)
{
  rf <- 0.5*(1-exp(-mar.sp/50))
  qtl.rf <- 0.5*(1-exp(-qtl.dist/50))
  n.chr <- length(n.mar)
  
  # possible errors in arguments
  if(length(mar.sp) != sum(n.mar) - n.chr) 
    stop("Length of mar.sp doesn't conform to n.mar.")
  if(length(qtl.chr) != length(qtl.mar) ||
     length(qtl.chr) != length(qtl.rf) ||
     length(qtl.chr) != length(qtl.eff)) 
    stop("Lengths of qtl.chr, qtl.mar, qtl.rf, qtl.eff must all be the same.")
  if(max(qtl.chr) > n.chr || min(qtl.chr) < 1) 
    stop("Entries qtl.chr must be between 1 and n.chr, inclusive.")

  # convert qtl.mar to cumulative numbers 0,1,2,...,sum(n.mar)-1
  qtl.mar <- cumsum(c(0,n.mar))[qtl.chr]+qtl.mar-1

  z <- .C("R_simbc_qtl",
          as.integer(n.ind),
          as.integer(n.chr),
          as.integer(n.mar),
          as.double(rf),
          geno=as.integer(rep(0,n.ind*sum(n.mar))),
          pheno=as.double(rep(0,n.ind)),
          as.integer(length(qtl.chr)),
          as.integer(qtl.chr),
          as.integer(qtl.mar),
          as.double(qtl.rf),
          as.double(qtl.eff),
          as.double(sigma),
          PACKAGE="qtlsim")

  list(geno=matrix(z$geno,ncol=n.chr*n.mar,nrow=n.ind),pheno=z$pheno)
}

# end of sims.R
