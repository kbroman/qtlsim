######################################################################
#
# util.R: Various utility functions
#
# copyright (c) 2001-6, Karl W Broman
# Oct, 2001; July, 2001, June, 2002; June, 2003; Oct 2006
# Licensed under the GNU General Public License version 3
#
# Part of the R/qtlsim package
#
# convertSim: convert output of simbc() to a cross object
#              as used by the R/qtl package
# which.pos
# which.correct
# which.correct.sub
# combine.sims
# convertOutput
#
######################################################################

convertSim <-
function(object, n.chr=9, chr.len=100)
{
  dat <- object
  n.mar <- ncol(dat$geno)/n.chr
  geno <- vector("list",n.chr)
  for(i in 1:n.chr) {
    geno[[i]]$data <- dat$geno[,n.mar*(i-1)+1:n.mar] + 1
    colnames(geno[[i]]$data) <-
      paste("D",i,"M",1:n.mar,sep="")
    geno[[i]]$map <- seq(0,chr.len,length=n.mar)
    names(geno[[i]]$map) <- colnames(geno[[i]]$data)
  }
  names(geno) <- as.character(1:n.chr)
  cross <- list(geno=geno,pheno=as.matrix(dat$pheno))
  colnames(cross$pheno) <- "pheno"
  class(cross) <- c("bc","cross")
  cross
}

which.pos <-
function(vect,cs)
{
  if(length(vect)==0) return(numeric(0))
  sapply(vect,function(a,b) {
    chr <- (1:length(b))[a <= b][1]
    c(chr,a-c(0,b)[chr])
  },cs)
}


which.correct.sub <-
function(infer,truth,within=1)
{
  if(length(infer)==0)
    return(rep(0,ncol(truth)+2))

  z <- .C("R_which_correct",
          as.integer(ncol(infer)),
          as.integer(infer[1,]),
          as.integer(infer[2,]),
          as.integer(ncol(truth)),
          as.integer(truth[1,]),
          as.integer(truth[2,]),
          as.integer(within),
          cor=as.integer(rep(0,ncol(truth))),
          incor=as.integer(c(0,0)),
          PACKAGE="qtlsim")
  c(z$cor,z$incor)
}



which.correct <-
function(results,
         truth=rbind(c(1,1,2,2,3,4,5),c(4,8,4,8,6,4,1)),
         within=1)
{
  flag <- 0
  if(!is.list(results[[1]])) {
    results <- list(results)
    flag <- 1
  }
  output <- lapply(results,function(a,b,d) {
      x <- t(sapply(a,which.correct.sub,b,d))
      n <- ncol(x)-2
      colnames(x) <- c(paste("QTL",1:n,sep=""),"incor.linked","incor.unlinked")
      x },truth,within)
  if(flag) output <- output[[1]]
  output
}


combine.sims <-
function(...)
{
  a <- list(...)
  if(length(a) < 2)
    stop("You must give at least two arguments.")
  if(length(unique(sapply(a,length))) != 1)
    stop("All arguments must be lists of the same length.")

  res <- a[[1]]
  for(i in 2:length(a))
    for(j in 1:length(res))
      res[[j]] <- c(res[[j]],a[[i]][[j]])

  res
}

combine.mcmc <-
function(...)
{
  a <- list(...)
  if(length(a) < 2)
    stop("You must give at least two arguments.")
  if(length(unique(sapply(a,length))) != 1)
    stop("All arguments must be lists of the same length.")

  res <- a[[1]]
  for(i in 2:length(a)) {
    res[[1]] <- combine.sims(res[[1]],a[[i]][[1]])
    res[[2]] <- rbind(res[[2]],a[[i]][[2]])
    res[[3]] <- rbind(res[[3]],a[[i]][[3]])
  }
  res
}

plot_res <-
function(x,...,ylim=NULL)
{
  tab <- sapply(x,apply,2,mean)
  par(...)
  if(is.null(ylim)) ylim <- c(0,1)
  if(ncol(tab) == 5) colors <- c("white","gray80","white","gray80","gray50")
  else colors <- c("white",rep("gray80",5),rep("gray50",3),"white")
  if(nrow(tab)==9) {
    ylab <- rep("Proportion observed",9)
    main <- c(paste("QTL",c("1a","1b","2a","2b","3","4","5")),
              "Linked extraneous","Unlinked extraneous")
  }
  else {
    ylab <- rep("",nrow(tab))
    main <- rep("",nrow(tab))
  }
  for(i in 1:nrow(tab)) {
    barplot(tab[i,],
            col=colors,ylim=ylim,main=main[i],ylab=ylab[i])
      abline(h=seq(0,1,by=0.1),lty=2)
  }

}


######################################################################
#
# convertOutput
#
# Convert the output from anal.leaps, anal.mcmc or perm
# to a form that may be used by the function which.correct():
# a matrix whose first row is the chromosome number and second
# row is the marker
#
######################################################################

convertOutput <-
function(object, n.mar=rep(11,9))
{
  output <- object
  if(is.list(output)) output <- output[[1]]
  if(length(output) == 0) return(numeric(0))
  output <- sort(output)

  if(any(output < 1 | output > sum(n.mar)))
    stop("Numbers is output are outside the range of the number of markers.")

  chr <- rep(1:length(n.mar),n.mar)[output]
  mar <- output-cumsum(c(0,n.mar))[chr]
  rbind(chr=chr,mar=mar)
}

# end of util.R
