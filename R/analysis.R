######################################################################
#
# analysis.R: R functions to perform QTL analysis
#
# copyright (c) 2001-3, Karl W Broman
# Nov, 2001; July, 2001; June 2003
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/qtlsim package
#
# anal.multi: Simulate backcross, perform ANOVA, CIM, and forward 
#             selection (with BIC and permutation tests) many times
# anal.multi2:Simulate backcross, perform ANOVA, CIM, forward 
#             selection (with BIC and permutation tests) and MCMC many times
# anal.all:   Perform ANOVA, CIM and forward selection on a backcross
# sim.null:   Simulate under the null hypothesis to get LOD thresholds
#             for ANOVA and CIM
# anal.leaps
# mcmc
# sim.mcmc
# 
######################################################################

anal.multi <-
function(n.sim=1000, cim.steps=c(3,5,7,9,11),bic.mult=c(2,2.5,3),
         max.steps=13, n.perm=1000, alpha=0.05, thresh, drop=1.5,
         n.ind=100, n.mar = rep(11,9), mar.sp = rep(10,10*9),
         qtl.chr = c(1,1,2,2,3,4,5), qtl.mar = c(4,8,4,8,6,4,1),
         qtl.dist = rep(0,7), qtl.eff = c(1,1,1,-1,1,1,1)/2,
         sigma=1)
{
  rf <- 0.5*(1-exp(-mar.sp/50))
  qtl.rf <- 0.5*(1-exp(-qtl.dist/50))
  n.chr <- length(n.mar)
  tot.mar <- sum(n.mar)
  n.cim <- length(cim.steps)
  cim.steps <- rev(sort(cim.steps))
  n.bic <- length(bic.mult)
  if(max(cim.steps) > max.steps) {
    warning("max.steps must be at least as big as max(cim.steps).")
    max.steps <- max(cim.steps)
  }
  # possible errors in arguments
  if(length(drop)==1) 
    drop <- rep(drop,n.cim+1)
  else
    if(length(drop) != n.cim+1)
      stop("Length of drop must be 1 or = 1+n.cim.")
  if(length(mar.sp) != sum(n.mar) - n.chr) 
    stop("Length of mar.sp doesn't conform to n.mar.")
  if(length(qtl.chr) != length(qtl.mar) ||
     length(qtl.chr) != length(qtl.rf) ||
     length(qtl.chr) != length(qtl.eff)) 
    stop("Lengths of qtl.chr, qtl.mar, qtl.rf, qtl.eff must all be the same.")
  if(max(qtl.chr) > n.chr || min(qtl.chr) < 1) 
    stop("Entries qtl.chr must be between 1 and n.chr, inclusive.")

  # thresholds
  if(!missing(thresh)) {
    if(length(thresh)==1) 
      thresh <- rep(thresh,n.cim+1)
    else
      if(length(thresh) != n.cim+1)
        stop("Length of thresh must be 1 or = 1+n.cim.")
  }
  else {
    if(length(n.mar) != 9 || any(n.mar != rep(11,9)) ||
       all(mar.sp != rep(10,10*9)) ||
       (n.ind != 100 && n.ind != 250 && n.ind != 500 && n.ind != 1000) ||
       any(is.na(match(cim.steps,c(3,5,7,9,11))))) {
      thresh <- rep(3,n.cim+1)
      warning("Using LOD threshold of 3.0")
    }
    else {
      if(n.ind == 100)
        thresh <- c(2.56,3.50,4.12,4.64,5.13,5.60)
      if(n.ind == 250) 
        thresh <- c(2.52,3.23,3.56,3.77,3.95,4.09)
      if(n.ind == 500)
        thresh <- c(2.50,3.15,3.38,3.51,3.60,3.67)
      if(n.ind == 1000)
        thresh <- c(2.48,3.08,3.27,3.38,3.43,3.47)
      thresh <- thresh[c(1,match(cim.steps,c(3,5,7,9,11))+1)]
    }
  }

  # convert qtl.mar to cumulative numbers 0,1,2,...,sum(n.mar)-1
  qtl.mar <- cumsum(c(0,n.mar))[qtl.chr]+qtl.mar-1

  z <- .C("R_anal_multi",
          as.integer(n.ind),
          as.integer(n.chr),
          as.integer(n.mar),
          as.double(rf),
          as.integer(n.sim),
          as.integer(length(qtl.chr)),
          as.integer(qtl.chr),
          as.integer(qtl.mar),
          as.double(qtl.rf),
          as.double(qtl.eff),
          as.double(sigma),
          as.integer(n.cim),
          as.integer(cim.steps),
          as.integer(max.steps),
          as.integer(n.bic),
          as.double(bic.mult),
          as.double(thresh),
          as.double(drop),
          as.integer(n.perm),
          as.integer(ceiling(alpha*n.perm)),
          n.qtl=as.integer(rep(0,n.sim*(2+n.cim+n.bic))),
          chr.id=as.integer(rep(0,n.sim*tot.mar*(2+n.cim+n.bic))),
          mar.id=as.integer(rep(0,n.sim*tot.mar*(2+n.cim+n.bic))),
          as.integer(rep(0,tot.mar*(n.ind+2))),
          as.double(rep(0,(tot.mar+2)^2+tot.mar+1+n.ind*(tot.mar+4))),
          PACKAGE="qtlsim")

  n.qtl <- matrix(z$n.qtl,ncol=2+n.cim+n.bic)
  id <- array(c(z$chr.id,z$mar.id),c(tot.mar,n.sim,2+n.cim+n.bic,2))
  results <- vector("list",2+n.cim+n.bic)
  for(i in 1:(2+n.cim+n.bic)) {
    results[[i]] <- vector("list",n.sim)
    for(j in 1:n.sim) {
      if(!n.qtl[j,i]) results[[i]][[j]] <- numeric(0)
      else {
        results[[i]][[j]] <- t(id[1:n.qtl[j,i],j,i,])
        if(n.qtl[j,i]==1)
          results[[i]][[j]] <- matrix(results[[i]][[j]],ncol=1)
        rownames(results[[i]][[j]]) <- c("chr","mar")
      }
    }
  }
  if(n.cim > 1) results[2:(n.cim+1)] <- results[(n.cim+1):2]
  names(results) <- c("ANOVA",paste("cim:",rev(cim.steps),sep=""),
                      paste("bic:",bic.mult,sep=""),"perm")

  results
}


anal.multi2 <-
function(n.sim=1000, cim.steps=c(3,5,7,9,11),bic.mult=c(2,2.5,3),
         max.steps=13, n.perm=1000, alpha=0.05, thresh, drop=1.5,
         n.ind=100, n.mar = rep(11,9), mar.sp = rep(10,10*9),
         n.mcmc=1000,mcmc.bic=2.56,
         qtl.chr = c(1,1,2,2,3,4,5), qtl.mar = c(4,8,4,8,6,4,1),
         qtl.dist = rep(0,7), qtl.eff = c(1,1,1,-1,1,1,1)/2,
         sigma=1)
{
  rf <- 0.5*(1-exp(-mar.sp/50))
  qtl.rf <- 0.5*(1-exp(-qtl.dist/50))
  n.chr <- length(n.mar)
  tot.mar <- sum(n.mar)
  n.cim <- length(cim.steps)
  cim.steps <- rev(sort(cim.steps))
  n.bic <- length(bic.mult)
  if(max(cim.steps) > max.steps) {
    warning("max.steps must be at least as big as max(cim.steps).")
    max.steps <- max(cim.steps)
  }
  # possible errors in arguments
  if(length(drop)==1) 
    drop <- rep(drop,n.cim+1)
  else
    if(length(drop) != n.cim+1)
      stop("Length of drop must be 1 or = 1+n.cim.")
  if(length(mar.sp) != sum(n.mar) - n.chr) 
    stop("Length of mar.sp doesn't conform to n.mar.")
  if(length(qtl.chr) != length(qtl.mar) ||
     length(qtl.chr) != length(qtl.rf) ||
     length(qtl.chr) != length(qtl.eff)) 
    stop("Lengths of qtl.chr, qtl.mar, qtl.rf, qtl.eff must all be the same.")
  if(max(qtl.chr) > n.chr || min(qtl.chr) < 1) 
    stop("Entries qtl.chr must be between 1 and n.chr, inclusive.")

  # thresholds
  if(!missing(thresh)) {
    if(length(thresh)==1) 
      thresh <- rep(thresh,n.cim+1)
    else
      if(length(thresh) != n.cim+1)
        stop("Length of thresh must be 1 or = 1+n.cim.")
  }
  else {
    if(length(n.mar) != 9 || any(n.mar != rep(11,9)) ||
       all(mar.sp != rep(10,10*9)) ||
       (n.ind != 100 && n.ind != 250 && n.ind != 500 && n.ind != 1000) ||
       any(is.na(match(cim.steps,c(3,5,7,9,11))))) {
      thresh <- rep(3,n.cim+1)
      warning("Using LOD threshold of 3.0")
    }
    else {
      if(n.ind == 100)
        thresh <- c(2.56,3.50,4.12,4.64,5.13,5.60)
      if(n.ind == 250) 
        thresh <- c(2.52,3.23,3.56,3.77,3.95,4.09)
      if(n.ind == 500)
        thresh <- c(2.50,3.15,3.38,3.51,3.60,3.67)
      if(n.ind == 1000)
        thresh <- c(2.48,3.08,3.27,3.38,3.43,3.47)
      thresh <- thresh[c(1,match(cim.steps,c(3,5,7,9,11))+1)]
    }
  }

  # convert qtl.mar to cumulative numbers 0,1,2,...,sum(n.mar)-1
  qtl.mar <- cumsum(c(0,n.mar))[qtl.chr]+qtl.mar-1

  z <- .C("R_anal_multi2",
          as.integer(n.ind),
          as.integer(n.chr),
          as.integer(n.mar),
          as.double(rf),
          as.integer(n.sim),
          as.integer(length(qtl.chr)),
          as.integer(qtl.chr),
          as.integer(qtl.mar),
          as.double(qtl.rf),
          as.double(qtl.eff),
          as.double(sigma),
          as.integer(n.cim),
          as.integer(cim.steps),
          as.integer(max.steps),
          as.integer(n.bic),
          as.double(bic.mult),
          as.double(thresh),
          as.double(drop),
          as.integer(n.perm),
          as.integer(ceiling(alpha*n.perm)),
          as.integer(n.mcmc),
          as.double(mcmc.bic*log(n.ind)),
          n.qtl=as.integer(rep(0,n.sim*(3+n.cim+n.bic))),
          chr.id=as.integer(rep(0,n.sim*tot.mar*(3+n.cim+n.bic))),
          mar.id=as.integer(rep(0,n.sim*tot.mar*(3+n.cim+n.bic))),
          as.integer(rep(0,tot.mar*(n.ind+4)+n.mcmc+4)),
          as.double(rep(0,(tot.mar+2)^2+n.ind*(tot.mar+4)+tot.mar+
                        n.mcmc+2)),
          PACKAGE="qtlsim")

  n.qtl <- matrix(z$n.qtl,ncol=3+n.cim+n.bic)
  if(any(n.qtl<0)) {
    warning("Some functions returned a negative number of QTLs identified\n")
  }
  id <- array(c(z$chr.id,z$mar.id),c(tot.mar,n.sim,3+n.cim+n.bic,2))
  results <- vector("list",3+n.cim+n.bic)
  for(i in 1:(3+n.cim+n.bic)) {
    results[[i]] <- vector("list",n.sim)
    for(j in 1:n.sim) {
      if(n.qtl[j,i]<1) results[[i]][[j]] <- numeric(0)
      else {
        results[[i]][[j]] <- t(id[1:n.qtl[j,i],j,i,])
        if(n.qtl[j,i]==1)
          results[[i]][[j]] <- matrix(results[[i]][[j]],ncol=1)
        rownames(results[[i]][[j]]) <- c("chr","mar")
      }
    }
  }
  if(n.cim > 1) results[2:(n.cim+1)] <- results[(n.cim+1):2]
  names(results) <- c("ANOVA",paste("cim:",rev(cim.steps),sep=""),
                      paste("bic:",bic.mult,sep=""),"perm","mcmc")

  results
}



anal.all <-
function(dat,cim.steps=7,max.steps=20)
{  
  gen <- dat$geno
  phe <- dat$pheno
  n.ind <- length(phe)
  tot.mar <- ncol(gen)
  if(max.steps > tot.mar) {
    warning("max.steps cannot be bigger than the number of markers.")
    max.steps <- tot.mar
  }
  if(cim.steps > max.steps) {
    warning("max.steps must be at least as big as cim.steps.")
    max.steps <- cim.steps
  }
    
  if(nrow(gen) != n.ind)
    stop("Number of rows in geno must equal the length of pheno.")
  
  z <- .C("R_anal_all",
          as.integer(n.ind),
          as.integer(tot.mar),
          as.integer(gen),
          as.double(phe),
          as.double(rep(0,(tot.mar+2)^2)), # workspace
          as.integer(cim.steps),
          as.integer(max.steps),
          lod=as.double(rep(0,tot.mar)),
          lodcim=as.double(rep(0,tot.mar)),
          index=as.integer(rep(0,tot.mar+1)),
          rss=as.double(rep(0,tot.mar+1)),
          PACKAGE="qtlsim")

  list(lod=z$lod,lodcim=z$lodcim,
       forw=cbind(index=c(0,z$index[1:max.steps]),rss=z$rss[1:(max.steps+1)]))
}


perm <-
function(dat, n.perm=1000, alpha=0.05)
{
  gen <- dat$geno
  phe <- dat$pheno
  n.ind <- length(phe)
  tot.mar <- ncol(gen)
  if(nrow(gen) != n.ind)
    stop("Number of rows in geno must equal the length of pheno.")
  
  z <- .C("R_forw_perm",
          as.integer(n.ind),
          as.integer(tot.mar),
          as.integer(gen),
          as.double(phe),
          index=as.integer(rep(0,tot.mar)),
          as.integer(n.perm),
          as.integer(ceiling(alpha*n.perm)),
          n.chosen=as.integer(0),
          as.double(rep(0,n.ind*(tot.mar+3))),
          PACKAGE="qtlsim")
  
  if(z$n.chosen==0) return(numeric(0))
  return(result=z$index[1:z$n.chosen])
}


sim.null <-
function(n.ind=100, n.mar=rep(11,9),mar.sp=rep(10,10*9),
         cim.steps=c(3,5,7,9,11),n.sim=10000)
{
  rf <- 0.5*(1-exp(-mar.sp/50))
  n.chr <- length(n.mar)
  tot.mar <- sum(n.mar)
  if(length(mar.sp) != tot.mar-n.chr)
    stop("Length of mar.sp doesn't conform to n.mar.")
  
  n.cim <- length(cim.steps)
  cim.steps <- rev(sort(cim.steps))

  z <- .C("R_sim_null",
          as.integer(n.ind),
          as.integer(n.chr),
          as.integer(n.mar),
          as.integer(tot.mar),
          as.double(rf),
          as.integer(n.cim),
          as.integer(cim.steps),
          as.integer(n.sim),
          maxlod=as.double(rep(0,n.sim*(n.cim+1))),
          as.integer(rep(0,tot.mar*(n.ind+1))),
          as.double(rep(0,n.ind+(tot.mar+2)^2+2*tot.mar+1)),
          PACKAGE="qtlsim")

  z <- matrix(z$maxlod,nrow=n.sim)
  colnames(z) <- c("ANOVA",paste("CIM:",cim.steps,sep=""))
  z[,c(1,seq(ncol(z),2,by=-1))]
}

#id <-
#function(lod, thresh=3, drop=2.2, n.mar=rep(11,9))
#{
#  n.chr <- length(n.mar)
#  tot.mar <- sum(n.mar)
#  if(length(lod) != tot.mar)
#    stop("Length of lod must equal sum(n.mar).")
#
#  z <- .C("R_identify_qtls",
#          as.integer(n.chr),
#          as.integer(n.mar),
#          as.double(lod),
#          as.double(thresh),
#          as.double(drop),
#          n.qtl=as.integer(0),
#          chr.id=as.integer(rep(0,tot.mar)),
#          mar.id=as.integer(rep(0,tot.mar)),
#          as.integer(rep(0,tot.mar)),
#          as.double(rep(0,tot.mar)),
#          PACKAGE="qtlsim")
#
#  if(!z$n.qtl) return(NULL)
#  rbind(chr=z$chr.id[1:z$n.qtl],mar=z$mar.id[1:z$n.qtl])
#
#}  



#id.bic <-
#function(forw.res, mult=2, n.mar=rep(11,9),n.ind=100)
#{
#  n.chr <- length(n.mar)
#  tot.mar <- sum(n.mar)
#  rss <- forw.res[,2]
#  index <- forw.res[-1,1]
#
#  z <- .C("R_identify_qtls_bic",
#          as.integer(n.chr),
#          as.integer(n.mar),
#          as.integer(n.ind),
#          as.integer(index),
#          as.double(rss),
#          as.integer(length(index)),
#          n.qtl=as.integer(0),
#          chr.id=as.integer(rep(0,length(index))),
#          mar.id=as.integer(rep(0,length(index))),
#          as.double(mult),
#          PACKAGE="qtlsim")
#          
#  if(!z$n.qtl) return(NULL)
#  rbind(chr=z$chr.id[1:z$n.qtl],mar=z$mar.id[1:z$n.qtl])
#}  


anal.leaps <-
function(dat, method=c("forward","backward","forwback",
                "forwardleap","backwardleap"),
         bic.mult=2.5, max.steps=30)
{
  method <- match.arg(method)

  x <- dat$geno
  y <- dat$pheno
  bic0 <- log(sum((y-mean(y))^2))
  n <- length(y)

  id <- 1:ncol(x)

  if(method == "forward" || method=="backward") {
    out <- regsubsets(x=x,y=y, method=method,nvmax=max.steps)
    out <- summary(out)
  }
  else if(method=="forwback") {
    out1 <- regsubsets(x=x,y=y,method="forward",nvmax=max.steps)
    out1 <- summary(out1)
    x <- x[,out1$which[max.steps,-1]]
    id2 <- id[out1$which[max.steps,-1]]
    out2 <- regsubsets(x=x,y=y,method="backward",nvmax=max.steps)
    out2 <- summary(out2)
  }
  else {
    method1 <- substr(method,0,nchar(method)-4)
    out <- regsubsets(x=x,y=y,method=method1,nvmax=max.steps)
    out <- summary(out)
    x <- x[,out$which[max.steps,-1]]
    id <- id[out$which[max.steps,-1]]
    out <- regsubsets(x=x,y=y,method="exhaustive",nvmax=max.steps)
    out <- summary(out)
  }
  if(method != "forwback") {
    bic <- log(out$rss)+(1:length(out$rss))*bic.mult*log(n)/n
    o <- (1:length(bic))[bic==min(bic)][1]
  }
  else { # for forw/backw method, look at *all* models fit
    bic1 <- log(out1$rss)+(1:length(out1$rss))*bic.mult*log(n)/n
    bic2 <- log(out2$rss)+(1:length(out2$rss))*bic.mult*log(n)/n
    if(min(bic1) < min(bic2)) {
      bic <- bic1
      out <- out1
    }
    else {
      bic <- bic2
      out <- out2
      id <- id2
    }
    o <- (1:length(bic))[bic==min(bic)][1]
  }

  if(bic0 < min(bic)) result <- list(numeric(0),bic0)
  else result <- list(id[out$which[o,-1]],min(bic))
  
  names(result) <- c("id","bic")
  result
}

anal.mcmc <-
function(dat,bic.mult=2.5,n.steps=1000,
          start=NULL)
{
  gen <- dat$geno
  phe <- dat$pheno
  n.ind <- nrow(gen)
  tot.mar <- ncol(gen)
  if(length(phe) != n.ind)
    stop("length(dat$pheno) must equal nrow(dat$geno).")

  if(is.null(start))
    indicate <- rep(0,tot.mar+1)
  else {
    indicate <- c(length(start),rep(0,tot.mar))
    indicate[start+1] <- 1
  }

  z <- .C("R_mcmc_ms",
          as.integer(n.ind),
          as.integer(tot.mar),
          as.integer(gen),
          as.double(phe),
          as.double(rep(0,(tot.mar+2)^2)), # xpx
          as.integer(n.steps),
          n.qtl.id=as.integer(0),
          qtl.id=as.integer(rep(0,tot.mar)),
          neg.log.post=as.double(0),
          first.seen=as.integer(0),
          as.integer(rep(0,tot.mar+1)),
          as.integer(indicate),
          as.double(bic.mult*log(n.ind)), # delta
          n.qtl.id.list=as.integer(rep(0,n.steps+1)),
          post.list=as.double(rep(0,n.steps+1)),
          PACKAGE="qtlsim")
  
  if(z$n.qtl.id==0) id <- numeric(0)
  else id <- z$qtl.id[1:z$n.qtl.id]
  list(id=id,neg.log.post=z$neg.log.post,first=z$first.seen,
       all.n.qtl=z$n.qtl.id.list,all.post=z$post.list)
}

          
    
                
sim.mcmc <-
function(n.sim=1000,max.steps=28,bic.mult=2.5, n.steps=1000,
         n.ind=100, n.mar = rep(11,9), mar.sp = rep(10,10*9),
         qtl.chr = c(1,1,2,2,3,4,5), qtl.mar = c(4,8,4,8,6,4,1),
         qtl.dist = rep(0,7), qtl.eff = c(1,1,1,-1,1,1,1)/2,
         sigma=1,only.two=FALSE)
{
  results <- vector("list",5)
  bic <- matrix(ncol=5,nrow=n.sim)
  mcmc.info <- matrix(nrow=n.sim,ncol=2)
  names(results) <- colnames(bic) <-
    c("forward","backward","forw.leap","back.leap","mcmc")
  colnames(mcmc.info) <- c("first","n.models")
  for(i in 1:5)
    results[[i]] <- vector("list",n.sim)
  cs <- cumsum(n.mar)
  for(i in 1:n.sim) {
    dat <- simbc(n.ind,n.mar,mar.sp,qtl.chr,qtl.mar,qtl.dist,
                  qtl.eff,sigma)
    if(only.two) {
      temp <- anal.leaps(dat,"forward",bic.mult,max.steps)
      results[[1]][[i]] <- which.pos(temp[[1]],cs)
      bic[i,1] <- temp[[2]]
    }
    else {
      for(j in 1:4) {
        temp <- anal.leaps(dat,c("forward","backward","forwardleap",
                                  "backwardleap")[j],bic.mult,max.steps)
        results[[j]][[i]] <- which.pos(temp[[1]],cs)
        bic[i,j] <- temp[[2]]
      }
    }
    temp <- anal.mcmc(dat,bic.mult,n.steps,NULL)
    results[[5]][[i]] <- which.pos(temp[[1]],cs)
    bic[i,5] <- temp[[2]]
    mcmc.info[i,1] <- temp[[3]]
    mcmc.info[i,2] <- length(unique(round(temp[[5]],9)))
  }
  if(only.two) return(list(results=results[c(1,5)],bic=bic[,c(1,5)],
                           mcmc.info=mcmc.info))
  list(results=results,bic=bic,mcmc.info=mcmc.info)
}

# end of analysis.R
