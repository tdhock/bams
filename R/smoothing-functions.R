fit.gada <- function
### Run the first fitting steps of the gada algorithm.
(pro
### Profile data.frame.
 ){
  require(gada)
  gen.info <- data.frame(probe=rownames(pro),
                         chr=pro$chromosome,
                         pos=pro$position)
  gada.obj <- setupGADAgeneral(pro$logratio,gen.info=gen.info)
  SBL(gada.obj,estim.sigma2=TRUE)
}

gada.results <- function
### Recover a matrix of smoothed signals from the gada results.
(pro,
### Profile data.frame.
 fit.list
### List of gada results. Each gada result is a list, one element for
### each chromosome.
 ){
  each.chrom(pro,function(chr){
    Y <- chr$logratio
    chrName <- as.character(chr$chromosome)[1]
    ## need to return mat[nparam,nprobes]
    do.call(rbind,lapply(fit.list,function(fit){
      names(fit) <- as.character(attr(fit,"chr"))
      seg.df <- WextIextToSegments(fit[[chrName]])
      for(i in 1:nrow(seg.df)){
        seg <- seg.df[i,]
        Y[seg$IniProbe:seg$EndProbe] <- seg$MeanAmp
      }
      Y
    }))
  })
}

run.pelt <- function
### Smooth a profile using the PELT algorithm.
(profile,
### A profile data.frame.
 penalty="SIC",
### character specifying the penalty to use.
 values=0,
### vector of penalty parameters to try.
 FUN=cpt.mean
### PELT function to use.
 ){
  require(changepoint)
  each.chrom(profile,function(chr){
    Y <- chr$logratio
    chr.mat <- do.call(rbind,lapply(values,function(value){
      ##cat(sprintf("profile %s chr %s parameter %s\n",
        ##          chr$profile.id[1],chr$chromosome[1],as.character(value)))
      tryCatch({
        result <- FUN(Y,method="PELT",penalty=penalty,pen.value=value)
        ends <- result@cpts
        starts <- c(1,ends[-length(ends)]+1)
        for(i in seq_along(ends)){
          end <- ends[i]
          start <- starts[i]
          Y[start:end] <- result@param.est$mean[i]
        }
        Y
      },error=function(e){
        rep(NA,length(Y))
      })
    }))
    to.replace <- is.na(chr.mat[,1])
    if(any(to.replace)){
      last.ok <- max(which(!to.replace))
      chr.mat[to.replace,] <-
        matrix(chr.mat[last.ok,],sum(to.replace),ncol(chr.mat),byrow=TRUE)
    }
    chr.mat
  })
### Matrix of smoothed profiles: nparam x nprobes.
}

dnacopy.smoothvec <- function
### Smooth a profile using DNAcopy.
(profile,
### A profile data.frame.
 var,
### Smoothness variable.
 vals,
### Smoothness values.
 ...
### Other arguments, passed to segment.
 ){
  require(DNAcopy)
  do.call(rbind,lapply(vals,function(v){
    set.seed(1)
    cna <- with(profile,CNA(logratio,chromosome,position))
    smoothed.cna <- smooth.CNA(cna)
    segmented.cna <- segment(smoothed.cna)
    args <- list(smoothed.cna,verbose=0,v,...)
    names(args)[3] <- var
    ##print(args)
    result <- do.call(segment,args)
    profile$smooth <- NA
    for(i in 1:nrow(result$output)){
      ri <- result$output[i,]
      probes.on.segment <- with(ri,{
        profile$chromosome==chrom &
        loc.start <= profile$position &
        loc.end   >= profile$position
      })
      profile$smooth[probes.on.segment] <- ri$seg.mean
    }
    profile$smooth
  }))
### Matrix of smoothed profiles: nparam x nprobes.
}

each.chrom <- function
### Apply a smoothing function independently to each chromosome of a
### profile.
(profile,
### Profile data.frame.
 FUN
### Function that will take a profile data.frame for one chromosome
### and return a smoothing matrix for that chromosome: nparam x
### nprobes.
 ){
  chroms <- split(profile,profile$chromosome,drop=TRUE)
  mats <- lapply(chroms,function(pro){
    if(any(diff(pro$position) <= 0))
      stop("profile positions must be in strictly increasing order")
    result <- FUN(pro)
    stopifnot(is.matrix(result))
    result
  })
  smooth.mat <- matrix(NA,nrow(mats[[1]]),nrow(profile))
  for(chr in names(mats)){
    smooth.mat[,profile$chromosome==chr] <- mats[[chr]]
  }
  smooth.mat
### Matrix of smoothed profiles for the entire profile: nparam x
### nprobes.
}

run.cghseg <- function
### Run cghseg maximum likelihood DP segmentation and pick the model
### using picker.
(profile,
### Profile data.frame.
 picker
### Function that gets arguments const.lines (a data.frame with one
### line for each segment in the maximum likelihood segmentation for a
### chrom), smoothed.mat (matrix of smoothed profile for a chrom: kmax
### x nprobes), Y (vector of logratio measurements for a chrom), kmax
### (maximum number of segments to consider), n (number of probes),
### and should return the chosen smoothing vector for the chrom.
 ){
  require(cghseg)
  each.chrom(profile,function(chr){
    Y <- chr$logratio
    pos.range <- range(chr$pos)
    n <- length(Y)
    ##if(chr$chromosome[1]=="4")browser()
    kmax <- min(20,n)#k is the number of SEGMENTS not BREAKPOINTS
    cat(sprintf("chrom=%s kmax=%d\n",chr$chromosome[1],kmax))
    result <- cghseg:::segmeanCO(Y,kmax)
    const.lines <- do.call(rbind,lapply(1:kmax,function(k){
      th <- result$t.est[k, 1:k]
      rupt <- matrix(ncol = 2, c(c(1, th[1:k - 1] + 1), th))
      d <- data.frame(k,begin = rupt[, 1], end = rupt[, 2],
                      mean = apply(rupt,1,function(z) mean(Y[z[1]:z[2]])))
      transform(d,
                position.begin=chr$position[begin],
                position.end=chr$position[end])
    }))
    ## matrix( smoothing levels x positions )
    eachk <- split(const.lines,const.lines$k)
    smoothed.mat <- do.call(rbind,lapply(eachk,function(d){
      for(i in 1:nrow(d))Y[d[i,"begin"]:d[i,"end"]] <- d[i,"mean"]
      Y
    }))
    best <- picker(smoothed.mat=smoothed.mat,
           const.lines=const.lines,
           Y=Y,n=n,kmax=kmax,l=pos.range[2]-pos.range[1])
    ##PROBLEM: what if chromosomes are discontiguous?? We NEED TO take
    ##the profile data frame apart more carefully.
    best
  })
### Smoothing matrix nparam x nprobes.
}

### Iterative isotonic regression.
iir.basic <- function(x,y,K){
  n <- length(x)
  res <- matrix(rep(0,n*K),ncol=K)
  u <- rep(0,n) ; b <- rep(0,n)
  for(k in 1:K){
    res.pava<-pava(y-b,long.out=T)
    u <- res.pava$y
    res.pava <- pava(y-u,long.out=T,decreasing=T)
    b <- res.pava$y
    res[,k] <- u+b
  }
return(res)
}

## run.iir <- function(pro,...){
## ### from
## ###http://www.sites.univ-rennes2.fr/laboratoire-statistique/JEGOU/iir_1.0.tar.gz
##   require(iir)
##   each.chrom(pro,function(chr){
##     if(chr$chromosome=="Y")return(t(chr$logratio))
##     res <- with(chr,iir::iir(position,logratio,...))
##     t(res$fitted)
##   })
## }

### Run Wild Binary Segmentation.
run.wbs <- function(pro, th.const.vals){
  require(wbs)
  each.chrom(pro, function(chr){
    Y <- chr$logratio
    mean.mat <- matrix(NA, length(th.const.vals), length(Y))
    if(nrow(chr) >= 2){
      set.seed(1)
      w <- wbs(Y)
      list.or.na <- changepoints(w, th.const=th.const.vals)$cpt.th
      if(is.list(list.or.na)){
        for(param.i in seq_along(list.or.na)){
          changes <- sort(list.or.na[[param.i]])
          ends <- c(changes, length(Y))
          starts <- c(1, changes+1)
          for(seg.i in seq_along(starts)){
            start <- starts[[seg.i]]
            end <- ends[[seg.i]]
            mean.mat[param.i, start:end] <- mean(Y[start:end])
          }
        }
      }
    }
    mean.mat
  })
}


### This is a list of functions, each of which must return a matrix of
### smoothed profiles. The first argument of each function is a
### data.frame that represents a copy number profile, with at least
### columns: position logratio chromosome. We assume that positions
### are already sorted in ascending order p_1 < p_2. The second
### argument to each of these functions should be a vector of
### smoothing parameters, and there should be a default value. The
### matrix returned has 1 row for each parameter, and 1 column for
### each position.
smoothers <-
  list(stat.ergo.alpha=function(profile,lambda){
    
  },flsa=function(profile,
l2vals=c(1e-5,1e-4,1e-3,5e-3,1e-2,2e-2,5e-2,
                 10^seq(-1,2,l=100),
                 2e2,5e2,1e3,1e4,1e5,
                 1e6,1e7,1e8,1e9,1e10,1e11,1e12)
         ){
    require(flsa)
    each.chrom(profile,function(d){
      model <- flsa(d$logratio)
      flsaGetSolution(model,lambda2=l2vals)
    })
  },iir.aicc=function(pro,x=1){
    run.iir(pro,criterion="aicc")
  },iir.aic=function(pro,x=1){
    run.iir(pro,criterion="aic")
  },iir.bic=function(pro,x=1){
    run.iir(pro,criterion="bic")
  },iir.gcv=function(pro,x=1){
    run.iir(pro,criterion="gcv")
  },pelt.default=function(pro,values=0){
    run.pelt(pro)
  },pelt.meanvar.n=function(pro,values=sprintf("n*%.10f",10^seq(-2,1,l=100))){
    run.pelt(pro,"Manual",values,cpt.meanvar)
  },pelt.meanvar=function(pro,values=10^seq(0,3,l=100)){
    run.pelt(pro,"Manual",values,cpt.meanvar)
  },pelt.n=function(pro,values=sprintf("n*%.10f",10^seq(-8,1,l=100))){
    run.pelt(pro,"Manual",values)
  },fpop=function(pro,lambda=10^seq(-8,1,l=100)){
    require(fpop)
    each.chrom(pro, function(chr){
      Y <- chr$logratio
      beta <- length(Y) * lambda
      chr.mat <- matrix(NA, length(beta), length(Y))
      for(beta.i in seq_along(beta)){
        res <- Fpop(chr$logratio, beta[[beta.i]])
        ends <- res$t.est
        starts <- c(1, 1+res$t.est[-length(res$t.est)])
        for(seg.i in seq_along(starts)){
          start <- starts[[seg.i]]
          end <- ends[[seg.i]]
          chr.mat[beta.i, start:end] <- mean(Y[start:end])
        }
      }
      chr.mat
    })
  },multiBinSeg=function(pro, lambda=10^c(seq(-8,1,l=100))){
    require(fpop)
    each.chrom(pro, function(chr){
      Y <- chr$logratio
      n <- length(Y)
      lambda.mat <- matrix(NA, length(lambda), n)
      if(n > 1){
        max.breaks <- min(20,n-1)#k is the number of BREAKPOINTS not segments
        result <- multiBinSeg(Y, max.breaks)
        model.mat <- matrix(NA, n, max.breaks+1)
        ## The first column is the constant model.
        model.mat[,1] <- mean(Y)
        breaks <- result$t.est
        for(model.i in seq_along(breaks)){
          changes <- sort(breaks[1:model.i])
          ends <- c(changes, length(Y))
          starts <- c(1, changes+1)
          for(seg.i in seq_along(starts)){
            start <- starts[[seg.i]]
            end <- ends[[seg.i]]
            model.mat[start:end, model.i+1] <- mean(Y[start:end])
          }
        }
        ## apply(model.mat, 2, function(yhat)sum(diff(yhat)!=0))
        residual.mat <- Y - model.mat
        mean.sq.residual <- colMeans(residual.mat * residual.mat)
        n.segments <- seq_along(mean.sq.residual)
        for(lambda.i in seq_along(lambda)){
          l <- lambda[[lambda.i]]
          crit <- n.segments * l + mean.sq.residual
          opt.segments <- which.min(crit)
          lambda.mat[lambda.i, ] <- model.mat[, opt.segments]
        }
      }
      lambda.mat
    })
  }, wbs.default=function(pro, th.const=1.3){
    run.wbs(pro, th.const)
  },wbs.th.const=function(pro, th.const=seq(0, 10, by=0.05)){
    run.wbs(pro, th.const)
  },pelt.diffparam=function(pro,
      values=sprintf("diffparam*%.10f",10^seq(-2,1,l=100))){
    run.pelt(pro,"Manual",values)
  },pelt.manual=function(pro,values=10^seq(-2,1,l=100)){
    run.pelt(pro,"Manual",values)
  },pelt.asymptotic=function(pro,
### for some profiles, for values <= 8e-2, we see

###Error in PELT.mean.norm(data, value) : 
###  NA/NaN/Inf in foreign function call (arg 4)
###In addition: Warning message:
###In log(log((1 - alpha + exp(-2 * (pi^(1/2)) * exp(blogn/alogn)))^(-1/(2 *  :
###   NaNs produced

### so lets check for this, and use the closest value if this happens.

      values=c(1-10^seq(-20,-1,l=80),8*(10^seq(-1,-3,l=20)))){
    run.pelt(pro,"Asymptotic",values)
  },flsa.norm=function(profile,
alpha2.vals=c(1e-5,1e-4,1e-3,5e-3,1e-2,2e-2,5e-2,
                 10^seq(-1,2,l=100),
                 2e2,5e2,1e3,1e4,1e5,
                 1e6,1e7,1e8,1e9,1e10,1e11,1e12)
         ){
    require(flsa)
    each.chrom(profile,function(d){
      model <- flsa(d$logratio)
      l <- max(d$position)-min(d$position)
      flsaGetSolution(model,lambda2=nrow(d)*1e6*alpha2.vals/l)
    })
  },cghseg.mBIC=function(profile,param=1){
    require(cghseg)
  each.chrom(profile,function(chr){
    Y <- chr$logratio
    if(length(Y)<10)return(rbind(rep(mean(Y),length(Y))))
    cghdata <- new("CGHdata",Y)
    cghoptions <- new("CGHoptions")
    result <- uniseg(cghdata,cghoptions)
    sm <- Y
    segs <- result@mu$Y
    for(i in 1:nrow(segs)){
      sm[segs[i,"begin"]:segs[i,"end"]] <- segs[i,"mean"]
    }
    rbind(sm)
  })
},cghseg.lambda=function(profile,lambda=10^c(seq(-8,-3,l=10),
                                  seq(-2,3,l=100),
                                  seq(4,8,l=10))){
## The point of this function is to try to recover the best
## segmentation according to some smart weighting of model complexity
## and signal reconstruction error. First, I wrote the following to
## Guillem
  
## Hi Guillem. I experimented a bit with the cghseg package this morning,
## and I would like some additional advice, please!

## I found the cghseg:::segmeanCO function which returns the t.est matrix
## of indices of the optimal segments for each k. This is useful, but I
## need a criterion besides k for selecting the level of smoothness.

## To provide some context, I have annotations about where the
## breakpoints really are, and I would like to find a smoothing function
## that agrees with those annotations. For some profiles there are no
## breakpoints, and for some profiles, there are many breakpoints. Thus
## fixing a k and applying that to each profile will certainly not agree
## with all profiles.

## In DNAcopy there is a parameter "undo.SD" which merges adjacent
## segments if the difference of their means is less than the value you
## specify. By varying this parameter I can find a sweet spot that
## minimizes false positive and false negative rates with respect to my
## annotations. I suppose I could do something similar with the output of
## segmeanCO(), but do you have any code already written that does this?

## In GLAD there is a parameter "lambdabreak" which essentially picks the
## number of breakpoints for a smoothing x by minimizing

## error(x) + lambdabreak * complexity(x)

## This is nice because I can again try various levels of lambdabreak and
## just pick the one that gives me a function that maximizes the
## agreement with my annotation data.

## I looked at your definition of BIC in the Kselection() function, and I
## was thinking of adapting it for this purpose, but I was unsure about
## which terms correspond to error and which correspond to model
## complexity. Can you clarify?

##         if (SSbg == 0) {
##             BIC[k] = gamma.coef + 0.5 * k * log(SSall) - 0.5 *
##                 sum(log(nk)) + (0.5 - k) * log(n)
##         }
##         else {
##             BIC[k] = 0.5 * (n - k + 1) * log(1 + (SSbg)/(SSall -
##                 SSbg)) + gamma.coef + 0.5 * k * log(SSall) -
##                 0.5 * sum(log(nk)) + (0.5 - k) * log(n)
##         }

######################### Then, he wrote back:
  
## I guess the simplest approach is to recover all best segmentations
## in 1 to K segments (with segmeanCO), then compute the agreement
## with your annotation data for all k and finally pick the k that
## maximizes the agreement with your annotation data.  That would be
## very similar to your lambda-break strategy with GLAD (k playing the
## role of the complexity and the likelihood or L2 loss playing the
## role of the error).

## I think, it might be more difficult to play with the mBIC criteria
## since the mBIC penalization takes into account many different type
## of complexity (number of parameters, length of segments, number of
## possible segmentations) In the computation of the BIC, (SSall -
## SSbg) is the L2 loss and can be interpreted as the error.  If you
## look at the case SSbg > 0, it could be rewritten as :

## BIC[k] =     0.5 * (n - k + 1) * ( log(SSall) - log(SSall - SSbg) ) 
##            + 0.5 * k * log(SSall) 
##            + gamma.coef
##            - 0.5 * sum(log(nk)) + (0.5 - k) * log(n)

##          = - 0.5 * (n - k + 1) * log(SSall - SSbg)
##            + 0.5 * (n + 1) *     log(SSall)
##            + gamma.coef
##            - 0.5 * sum(log(nk)) + (0.5 - k) * log(n)

## where
  
## a) gamma.coef penalize for the number of possible segmentation (n-1
## \choose k-1) ;
  
## b) 0.5 * sum(log(nk)) is the entropy of the segment lengths
## (segmentations with very smaller segments get higher penalty) ;
  
## c) (0.5 - k) * log(n) penalize the number of parameters to estimate
  
## d) finally the likelihood is multiplied by 0.5 * (n - k + 1).
  
## Due to this last term it is not easy to get a nice additive formula
## of the form : error(x) + complexity(x).

## (The case SSbg == 0 correspond to  log(SSall) - log(SSall - SSbg) == 0)

## Hopefully this will answer your questions. Let me know if that is
## not the case or if you have additional questions.  Cheers, Guillem
  run.cghseg(profile,function(kmax,Y,const.lines,n,smoothed.mat,...){
    lterms <- rep(NA,kmax)
    rterms <- rep(NA,kmax)
    ybar <- mean(Y,na.rm=TRUE)
    SSall <- sum((Y-ybar)^2,na.rm=TRUE)
    for(k in 1:kmax){
      gamma.coef <- lgamma(0.5*(n-k-1)) - lgamma(0.5*(n+1))
      mu <- const.lines[const.lines$k==k,]
      nk <- with(mu,end-begin+1) ### !!! ok
      SSbg <- sum(nk*(mu$mean-ybar)^2,na.rm=TRUE)
      rterms[k] <-
        gamma.coef + 0.5*k*log(SSall)-0.5 * sum(log(nk))+(0.5-k)*log(n)
      lterms[k] <- if (SSbg ==0){
        0
      } else {
        0.5*(n-k+1) * log(1+ (SSbg)/(SSall-SSbg))
      }
    } # end k
    do.call(rbind,lapply(lambda,function(l){
      smoothed.mat[which.min(lterms + l*rterms),]
    }))
  })
},cghseg.k=function(profile,lambda=10^c(seq(-8,1,l=100))){
### From Lavielle 2005 p.8, where he suggests to use average squared
### residual + beta*(number of segments)
  run.cghseg(profile,function(kmax,Y,smoothed.mat,...){
    average.residual <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    for(k in 1:kmax){
      average.residual[k] <- mean((Y-smoothed.mat[k,])^2)
      penalty[k] <- k
    }
    do.call(rbind,lapply(lambda,function(l){
      smoothed.mat[which.min(average.residual + l*penalty),]
    }))
  })
},cghseg.k.sqrt=function(profile,lambda=10^c(seq(-1,5,l=100))){
### Simulations indicate that a factor of sqrt(d) is needed to correct
### for number of points sampled, and 1/sqrt(l) is needed to correct
### for signal size in base pairs.
  run.cghseg(profile,function(kmax,Y,smoothed.mat,l,...){
    total.residual <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    d <- length(Y)
    for(k in 1:kmax){
      total.residual[k] <- sum((Y-smoothed.mat[k,])^2)
      penalty[k] <- k*sqrt(d/l)
    }
    do.call(rbind,lapply(lambda,function(l){
      smoothed.mat[which.min(total.residual + l*penalty),]
    }))
  })
},cghseg.k.var=function(profile,lambda=10^c(seq(-3,1,l=100))){
### Simulations indicate that a factor of s.hat^2 is necessary to
### correct for differences in scale of the samples.
  run.cghseg(profile,function(kmax,Y,smoothed.mat,...){
    total.residual <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    d <- length(Y)
    s.hat <- median(abs(diff(Y)))
    for(k in 1:kmax){
      total.residual[k] <- sum((Y-smoothed.mat[k,])^2)
      penalty[k] <- k*d*(s.hat^2)
    }
    do.call(rbind,lapply(lambda,function(l){
      smoothed.mat[which.min(total.residual + l*penalty),]
    }))
  })
},cghseg.k.sqrt.var=function(profile,lambda=10^c(seq(2,7,l=100))){
### Try using all 3 corrections shown to be beneficial in simulations:
### density d, base pairs length l, and scale s.hat.
  run.cghseg(profile,function(kmax,Y,smoothed.mat,l,...){
    total.residual <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    d <- length(Y)
    s.hat <- median(abs(diff(Y)))
    for(k in 1:kmax){
      total.residual[k] <- sum((Y-smoothed.mat[k,])^2)
      penalty[k] <- k*sqrt(d/l)*(s.hat^2)
    }
    do.call(rbind,lapply(lambda,function(l){
      smoothed.mat[which.min(total.residual + l*penalty),]
    }))
  })
},cghseg.k.sqrt.d=function(profile,lambda=10^c(seq(-3,0,l=100))){
### Use sqrt(d) instead of d in cghseg.k
  run.cghseg(profile,function(kmax,Y,smoothed.mat,...){
    total.residual <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    d <- length(Y)
    s.hat <- median(abs(diff(Y)))
    for(k in 1:kmax){
      total.residual[k] <- sum((Y-smoothed.mat[k,])^2)
      penalty[k] <- k*sqrt(d)
    }
    do.call(rbind,lapply(lambda,function(l){
      smoothed.mat[which.min(total.residual + l*penalty),]
    }))
  })
},cghseg.k.sqrt.d.var=function(profile,lambda=10^c(seq(-1,3,l=100))){
### Use the sqrt(d) correction with the scale correction.
  run.cghseg(profile,function(kmax,Y,smoothed.mat,...){
    total.residual <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    d <- length(Y)
    s.hat <- median(abs(diff(Y)))
    for(k in 1:kmax){
      total.residual[k] <- sum((Y-smoothed.mat[k,])^2)
      penalty[k] <- k*sqrt(d)*(s.hat^2)
    }
    do.call(rbind,lapply(lambda,function(l){
      smoothed.mat[which.min(total.residual + l*penalty),]
    }))
  })
},cghseg.k1=function(profile,lambda=10^c(seq(-8,1,l=100))){
### From Lavielle 2005 p.8, where he suggests to use average squared
### residual + beta*(number of segments)
  run.cghseg(profile,function(kmax,Y,smoothed.mat,...){
    average.residual <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    for(k in 1:kmax){
      average.residual[k] <- mean((Y-smoothed.mat[k,])^2)
      penalty[k] <- k - 1
    }
    do.call(rbind,lapply(lambda,function(l){
      smoothed.mat[which.min(average.residual + l*penalty),]
    }))
  })
},cghseg.k.less=function(profile,lambda=10^c(seq(-4,-1,l=10))){
### From Lavielle 2005 p.8, where he suggests to use average squared
### residual + beta*(number of segments)
  run.cghseg(profile,function(kmax,Y,smoothed.mat,...){
    average.residual <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    for(k in 1:kmax){
      average.residual[k] <- mean((Y-smoothed.mat[k,])^2)
      penalty[k] <- k
    }
    do.call(rbind,lapply(lambda,function(l){
      smoothed.mat[which.min(average.residual + l*penalty),]
    }))
  })
},cghseg.resid.corr=function(profile,lambda=10^c(seq(-8,1,l=100))){
### Use correletion coefficient of neighboring residuals instead of
### average squared residual.
  run.cghseg(profile,function(kmax,Y,smoothed.mat,...){
    resid.corr <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    for(k in 1:kmax){
      residuals <- Y-smoothed.mat[k,]
      resid.corr[k] <- cor(residuals[-1],residuals[-length(residuals)])
      penalty[k] <- k
    }
    resid.corr[is.na(resid.corr)] <- 0
    ##plot(penalty,resid.corr,type="l",ylim=c(-1,1));browser()
    do.call(rbind,lapply(lambda,function(l){
      smoothed.mat[which.min(resid.corr + l*penalty),]
    }))
  })
},cghseg.klogn=function(profile,lambda=10^c(seq(-8,1,l=100))){
### From Lavielle 2005 p.8, where he says the schwarz criterion means
### min average squared residual + 2 \sigma^2(\log n)/n * k
  run.cghseg(profile,function(kmax,Y,n,smoothed.mat,...){
    average.residual <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    for(k in 1:kmax){
      average.residual[k] <- mean((Y-smoothed.mat[k,])^2)
      penalty[k] <- k
    }
    do.call(rbind,lapply(lambda,function(l){
      smoothed.mat[which.min(average.residual + l*penalty*log(n)),]
    }))
  })
},cghseg.c=function(profile,cvals=10^c(seq(-8,-3,l=10),
                                  seq(-2,3,l=100),
                                  seq(4,8,l=10))){
### From Lavielle 2005 equation 14
  run.cghseg(profile,function(kmax,Y,n,smoothed.mat,...){
    average.residual <- rep(NA,kmax)
    kseq <- 1:kmax
    for(k in 1:kmax){
      average.residual[k] <- mean((Y-smoothed.mat[k,])^2)
    }
    s <- median(abs(diff(Y)))
    do.call(rbind,lapply(cvals,function(cval){
      pen <- 1+cval*log(n/kseq)
      smoothed.mat[which.min(average.residual + 2*s^2/n*kseq*pen),]
    }))
  })
},cghseg.schwarz=function(profile,sigma=10^c(seq(-8,-3,l=10),
                                  seq(-2,3,l=100),
                                  seq(4,8,l=10))){
### From Lavielle 2005 p.8, where he says the schwarz criterion means
### min average squared residual + 2 \sigma^2(\log n)/n * k, and lets
### learn the sigma parameter.
  run.cghseg(profile,function(kmax,Y,n,smoothed.mat,...){
    average.residual <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    for(k in 1:kmax){
      average.residual[k] <- mean((Y-smoothed.mat[k,])^2)
      penalty[k] <- k
    }
    do.call(rbind,lapply(sigma,function(s){
      smoothed.mat[which.min(average.residual + 2*s^2*log(n)/n*penalty),]
    }))
  })
},cghseg.schwarz.auto=function(profile,sigma=1){
### From Lavielle 2005 p.8, where he says the schwarz criterion means
### min average squared residual + 2 \sigma^2(\log n)/n * k. Here we
### estimate sigma using the median of differences.
  run.cghseg(profile,function(kmax,Y,n,smoothed.mat,...){
    average.residual <- rep(NA,kmax)
    penalty <- rep(NA,kmax)
    for(k in 1:kmax){
      average.residual[k] <- mean((Y-smoothed.mat[k,])^2)
      penalty[k] <- k
    }
    s <- median(abs(diff(Y)))
    stat <- average.residual + 2*s^2*log(n)/n*penalty
    smoothed.mat[which.min(stat),,drop=FALSE]
  })
},cghseg.min.abs.diff=function(profile,madvals=10^seq(-2,1,l=40)){
### The idea here is to do something like DNAcopy prune sd=...
  run.cghseg(profile,function(const.lines,smoothed.mat,...){
    eachk <- split(const.lines,const.lines$k)
    min.abs.diff <-
      sapply(eachk,with,min(abs(if(length(mean)==1)Inf else diff(mean))))
    do.call(rbind,lapply(madvals,function(MAD){
      smoothed.mat[which.min(min.abs.diff[min.abs.diff>MAD]),]
    }))
  })
},cghFLasso=function(profile,param=1){
  require(cghFLasso)
  t(cghFLasso(profile$logratio,profile$chromosome)$Esti.CopyN)
},glad.default=function(profile,param=1){
  runglad(profile,lambdabreak=8)
},glad.haarseg=function(profile,qvals=10^seq(-80,-1,l=30)){
  runglad(profile,breaksFdrQ=qvals,smoothfunc="haarseg")
},glad.lambdabreak=function(profile,
    lambdavals=c(1e-2,1,5,10,20,25,30,35,40,100,500,1e4,1e8)){
  runglad(profile,lambdabreak=lambdavals)
},glad.MinBkpWeight=function(profile,wvals=2^seq(-12,3,l=30)){
  runglad(profile,MinBkpWeight=wvals)
},gada=function(pro,Tvals=10^seq(-1,2,l=100)){
  fit <- fit.gada(pro)
  ## do all the elimination procedures in 1 go.
  be.list <- lapply(Tvals,function(Tval){
    BackwardElimination(fit,Tval,1)
  })
  gada.results(pro,be.list)
},gada.default=function(pro,param=0){
  fit <- fit.gada(pro)
  gada.results(pro,list(fit))
},dnacopy.default=function(profile,param=1){
  dnacopy.smoothvec(profile,"undo.SD",3)
},dnacopy.sd=function(profile,
    sdvals=c(seq(20,4,l=20),2^seq(2.5,-5,l=10)[-1],0)){
  dnacopy.smoothvec(profile,"undo.SD",sdvals,undo.splits="sdundo")
},dnacopy.prune=function(profile,prunevals=10^seq(-2,1,l=20)){
  dnacopy.smoothvec(profile,"undo.prune",prunevals,undo.splits="prune")
},dnacopy.alpha=function(profile,alphavals=10^seq(-80,-2,l=20)){
  dnacopy.smoothvec(profile,"alpha",alphavals)
})

runglad <- function
### Run glad to smooth a profile.
(profile,
### Profile data.frame.
 ...
### Smoothing parameter for glad.
 ){
  require(GLAD)
  L <- list(...)
  vals <- L[[1]]
  x <- transform(profile,
                 PosOrder=seq_along(position),
                 LogRatio=logratio,
                 PosBase=position,
                 Chromosome=chromosome)
  cgh <- as.profileCGH(x)
  smooth.vecs <- lapply(vals,function(v){
    args <- list(profileCGH=cgh,v)
    names(args)[2] <- names(L)[1]
    args <- c(args,L[-1])
    smooth <- tryCatch({
      args$OnlySmoothing <- TRUE
      do.call(daglad,args)
    },error=function(e){
      do.call(glad,args)
    })
    smooth$profileValues$Smoothing
  })
  do.call(rbind,smooth.vecs)
### Smoothing matrix nparam x nprobes.
}

several.breakpoints <- function
### Zero-one loss when breakpoint annotations are considered to mean
### presence of 1 or more breakpoints.
(counts,
### Counts of breakpoints of model.
 anns
### Annotations.
 ){
  ifelse(anns=="normal",
         ifelse(counts==0,0L,1L), ## normal
         ifelse(counts>0,0L,1L)) ## breakpoint
### Matrix of 0 and 1.
}

run.smoothers <- function
### Run several smoothers on a profile, quantifying agreement of the
### smoother with breakpoint annotations.
(profile,
### Profile data.frame.
 breakpoint.labels,
### Annotation data.frame.
 smooth.funs=smoothers,
### List of smoothing functions to apply to the profile.
 loss=several.breakpoints,
### Loss function of (breakpoint counts,labels) --- both matrices
### nparam x nlabels, and should return numeric matrix of same size.
 tosave=c("errors","seconds","parameters","breakpoint.labels"),
### Variables to save to the db directory.
 db=file.path(Sys.getenv("HOME"),"smooth")
### Location to save gzipped result files.
 ){
  stopifnot(is.data.frame(profile))
  counter <- function(x)sapply(x,function(i)sum(i==x))
  for(chrom in split(profile,profile$chromosome)){
    counts <- counter(chrom$position)
    repeats <- chrom[counts>1,]
    if(nrow(repeats)){
      print(repeats)
      stop("some repeated positions were found! most algorithms will break!")
    }
  }
  stopifnot(is.data.frame(breakpoint.labels))
  breakpoint.labels <- subset(breakpoint.labels,!is.na(annotation))
  if(nrow(breakpoint.labels)==0)stop("need at least 1 annotation")
  stopifnot(is.list(smooth.funs))
  if(is.null(names(smooth.funs)))stop("smooth.funs need to be named")
  if(any(nchar(names(smooth.funs))==0))stop("all smooth.funs need names")
  if(!file.exists(db))dir.create(db)
  smoothing.parameters <- lapply(smooth.funs,function(f)eval(formals(f)[[2]]))
  for(algorithm in names(smooth.funs)){
    print(algorithm)
    fun <- smooth.funs[[algorithm]]
    parameters <- smoothing.parameters[[algorithm]]
    seconds <- system.time({
      smoothed.profiles <- fun(profile)
    })["elapsed"]
    ## reality checks for result of smoothing function
    stopifnot(ncol(smoothed.profiles)==nrow(profile))
    stopifnot(nrow(smoothed.profiles)==length(parameters))

    bpts <- matrix(NA,length(parameters),nrow(breakpoint.labels))
    blabs <- matrix(NA,length(parameters),nrow(breakpoint.labels))
    pos <- profile$position
    for(i in 1:nrow(breakpoint.labels)){
      l <- breakpoint.labels[i,]
      probes <- l$min < pos & pos < l$max & profile$chromosome==l$chromosome
      smoothed.region <- smoothed.profiles[,probes,drop=FALSE]
      bpts[,i] <- apply(smoothed.region,1,function(x)sum(diff(x)!=0))
      blabs[,i] <- as.character(l$annotation)
    }
    print(bpts)

    errors <- loss(bpts,blabs)

    for(varname in tosave){
      if(varname %in% ls()){
        var <- get(varname)
        edir <- file.path(db,profile$profile.id[1],algorithm)
        if(!file.exists(edir))dir.create(edir,recursive=TRUE)
        this.base <- sprintf("%s.csv.gz",varname)
        efn <- file.path(edir,this.base)
        conn <- gzfile(efn,"w")
        write.table(var,conn,
                    row.names=FALSE,col.names=FALSE,quote=FALSE,sep=",")
        close(conn)
      }
    }
  }
### Nothing, the results are saved to files.
}

### Smoothing functions used in the article.
article.smoothers <-
  smoothers[c("cghseg.k","flsa.norm","flsa","glad.lambdabreak",
           "glad.MinBkpWeight","glad.haarseg",
           "cghseg.mBIC","cghFLasso","glad.default","dnacopy.default",
              "dnacopy.sd",
              "dnacopy.alpha",
              "dnacopy.prune"
              )]
