unseen.profile.error <- function
### Do n/t-fold cross-validation to estimate the error of global
### models with a small training set.
(stats,
### list with arrays for errors, false.positive, and false.negative
 prof.per.train
### t = approximate number of annotated profiles per training set.
 ){
  stat.names <- c("errors","false.positive","false.negative")
  nparam <- dim(stats$errors)[1]
  nprof <- dim(stats$errors)[2]
  narms <- dim(stats$errors)[3]
  nfolds <- floor(nprof/prof.per.train)
  test.error <- matrix(NA,nfolds,length(stat.names))
  colnames(test.error) <- stat.names
  fold <- sample(rep(1:nfolds,l=nprof))
  for(f in 1:nfolds){
    test.fold <- fold != f
    train.fold <- !test.fold
    train.err <- rep(NA,nparam) ## global model
    for(j in 1:nparam){
      train.err[j] <- mean(stats$errors[j,train.fold,],na.rm=TRUE)
    }
    global.picked <- pick.best.index(train.err)
    for(sn in stat.names){
      test.error[f,sn] <-
        mean(stats[[sn]][global.picked,test.fold,],na.rm=TRUE)
    }
    num.normal <- sum(stats$normal.anns[test.fold,])
    num.breakpoint <- sum(stats$breakpoint.anns[test.fold,])
    num.total <- num.normal+num.breakpoint
    FP <- test.error[f,"false.positive"]*num.total
    if(round(FP)-FP>1e-5){
      stop("FP not integral!")
    }
    test.error[f,"false.positive"] <- FP/num.normal
    test.error[f,"false.negative"] <-
      test.error[f,"false.negative"]*num.total/num.breakpoint
  }
  test.error
### matrix of estimated test errors, nfolds x 3
}

geom_tallrect <- function
### ggplot2 geom with xmin and xmax aesthetics that covers the entire
### y range.
(mapping=NULL,
 data=NULL,
 stat="identity",
 position="identity",
 ...){
  require(proto)
  require(grid)
  GeomTallRect <- proto(ggplot2:::GeomRect,{
    required_aes <- c("xmin", "xmax")
    draw <- draw_groups <- function(.,data,scales,coordinates,
                                    ymin=0,ymax=1,...){
      ymin <- unit(ymin,"npc")
      ymax <- unit(ymax,"npc")
      with(coord_transform(coordinates, data, scales),ggname(.$my_name(), {
        rectGrob(xmin, ymin, xmax - xmin, ymax-ymin,
                 default.units = "native", just = c("left", "bottom"), 
                 gp=gpar(
                   col=colour, fill=alpha(fill, alpha), 
                   lwd=size * .pt, lty=linetype, lineend="butt"
                   )
                 )
      }))
    }
  })
  GeomTallRect$new(mapping = mapping, data = data, stat = stat,
                   position = position, ...)
}

pick.best.index <- structure(function
### Minimizer for local models, described in article section 2.3
### "Picking the optimal model"
(err
### Vector of errors to minimize.
 ){
  nparam <- length(err)
  candidates <- which(err==min(err))
  if(length(err)==1)return(candidates)
  st <- abs(median(candidates)-candidates)
  middle <- candidates[which.min(st)]
  if(all(diff(err)==0))return(middle)
  if(nparam %in% candidates && 1 %in% candidates){
    cat("Warning: strange error profile, picking something near the center\n")
    print(as.numeric(err))
    d <- diff(candidates)>1
    if(any(d)){
      which(d)[1]
    }else{
      middle
    }
  }else if(1 %in% candidates){
    max(candidates)
  }else if(nparam %in% candidates){
    min(candidates)
  }else {
    middle
  }
### Integer index of the minimal error.
},ex=function(){
  stopifnot(pick.best.index(rep(0,100))==50)

  err <- rep(1,100)
  err[5] <- 0
  stopifnot(pick.best.index(err)==5)

  ## should pick the middle
  err <- rep(1,100)
  err[40:60] <- 0
  stopifnot(pick.best.index(err)==50)

  ## should pick the biggest
  err <- rep(1,100)
  err[1:60] <- 0
  stopifnot(pick.best.index(err)==60)

  ## should pick the smallest
  err <- rep(1,100)
  err[50:100] <- 0
  stopifnot(pick.best.index(err)==50)
})

estimate.test.error <- function
### Do leave-one-out cross-validation on chromosome arms.
(stats
### Named list with arrays errors, false.positive, false.negative,
### each of dim nparam x nprof x nfolds.
 ){
  stats <- stats[c("errors","false.positive","false.negative")]
  nparam <- dim(stats$errors)[1]
  nprof <- dim(stats$errors)[2]
  nfolds <- dim(stats$errors)[3]
  stat.names <- names(stats)
  local.loo <- matrix(NA,length(stats),nfolds)
  hybrid.loo <- matrix(NA,length(stats),nfolds)
  global.loo <- matrix(NA,length(stats),nfolds)
  rownames(global.loo) <- stat.names
  rownames(hybrid.loo) <- stat.names
  rownames(local.loo) <- stat.names
  train.err.mat <-
    matrix(NA,3,nfolds,dimnames=list(method=c("global","hybrid","local"),fold=NULL))
  for(fold in 1:nfolds){

    train.err <- rep(NA,nparam) ## global model
    for(j in 1:nparam){
      train.err[j] <- sum(stats$errors[j,,-fold],na.rm=TRUE)
    }
    ## save for hybrid approach
    global.train.err <-
      data.frame(train.err,param=1:length(train.err))
    global.picked <- pick.best.index(train.err)
    for(sn in stat.names){
      global.loo[sn,fold] <-
        mean(stats[[sn]][global.picked,,fold],na.rm=TRUE)
    }
    train.err.mat["global",fold] <- train.err[global.picked] ## for comparing train err

    ind.stats <- matrix(NA,nprof,length(stats))
    colnames(ind.stats) <- stat.names
    hybrid.train.errors <- rep(NA,nprof)
    for(i in 1:nprof){ ## hybrid models
      train.err <-
        apply(stats$errors[,i,-fold,drop=FALSE],1,sum,na.rm=TRUE)
      is.min <- train.err == min(train.err)
      global.subset <- global.train.err[is.min,]
      hybrid.picked <-
        with(global.subset,param[which.min(train.err)])
      hybrid.train.errors[i] <- train.err[hybrid.picked]
      for(sn in stat.names){ ## store test err for picked model
        ind.stats[i,sn] <- stats[[sn]][hybrid.picked,i,fold]
      }
    }
    hybrid.loo[,fold] <- colMeans(ind.stats,na.rm=TRUE)
    train.err.mat["hybrid",fold] <- sum(hybrid.train.errors)
    
    ind.stats <- matrix(NA,nprof,length(stats))
    colnames(ind.stats) <- stat.names
    local.train.errors <- rep(NA,nprof)
    for(i in 1:nprof){ ## local models
      train.err <-
        apply(stats$errors[,i,-fold,drop=FALSE],1,sum,na.rm=TRUE)
      local.picked <- pick.best.index(train.err)
      for(sn in stat.names){ ## store test err for picked model
        ind.stats[i,sn] <- stats[[sn]][local.picked,i,fold]
      }
      local.train.errors[i] <- train.err[local.picked]
    }
    local.loo[,fold] <- colMeans(ind.stats,na.rm=TRUE)
    train.err.mat["local",fold] <- sum(local.train.errors)
    
  }
  list(local=local.loo,
       hybrid=hybrid.loo,
       global=global.loo,
       train.err.mat=train.err.mat)
### Named list with elements local, hybrid, global, each a 3 x nfolds
### matrix.
}
