group.sofar <- function(x,rank=NULL,th=0.01){
  ## group the v matrix of so far fit according to
  ## the signs
  if(is.null(rank)){
    rank <- x$rank
  }
  # first few ranks to use
  stopifnot(rank<=ncol(x$V))
  g0 <- vector("list",rank)
  for(i in 1:length(g0)){
    v <- abs(x$V[,i])
    g <- rep(NA,length(v))
    g[abs(v)<th] <- 0
    g[abs(v)>=th] <- sign(x$V[abs(v)>=th,i])
    g0[[i]] <- g
  }
  interaction(g0)
}
fitted.sofar <- function(x,rank=NULL,center=NULL,log=F,scale=NULL,base=10){
  ## center: center to append
  ## log: T/F whether exponent
  ## scale: vector to scale the results
  ## base: log base for transformation
  y0 <- fitted.sofar.Z(x=x,rank=rank)
  
  y1 <- y0
  if(!is.null(scale)){
    stopifnot(ncol(y0)==length(scale))
    stopifnot(all(!is.na(scale)))
    y1 <- sweep(y0,2,scale,FUN="*")
  }
  
  y2 <- y1
  if(!is.null(center)){
    stopifnot(ncol(y1)==length(center))
    stopifnot(all(!is.na(center)))
    y2 <- sweep(y1,2,center,FUN="+")
  }
  
  y3 <- y2
  if(log){
    #zero <- y2==0
    y3 <- base^y2
    #y3[zero] <- 0
  }
  y3
}
fitted.sofar.Z.rank <- function(x,rank){
  ## derive fitted value on transformed scale by rank
  if(rank>x$rank){
    warning("rank ",rank," reset to ",x$rank,"\n")
    k <- x$rank
  }else{
    k <- rank
  }
  
  k <- as.integer(k)
  stopifnot(k>0)
  yhat <- with(x, U[,k,drop=F]%*%matrix(D[k],1,1)%*%t(V[,k,drop=F]))
  colnames(yhat) <- colnames(x$Y)
  v <-   x$V[,k,drop=T]
  yhat[,order(v)]
}
fitted.sofar.Z <- function(x,rank=NULL){
  ## derive fitted value on transformed scale
  if(is.null(rank)) {
    k <- x$rank
  }else{
    if(rank>x$rank){
      warning("rank ",rank," reset to ",x$rank,"\n")
      k <- x$rank
    }else{
      k <- rank
    }
  }
  k <- as.integer(k)
  stopifnot(k>0)
  sel <- seq(1,k)
  if(k==1){
    yhat <- with(x, U[,1,drop=F]%*%matrix(D[1],1,1)%*%t(V[,1,drop=F]))
  }else{
    yhat <- with(x, U[,sel]%*%diag(D[sel])%*%t(V[,sel]))
  }
  yhat
}
residuals.sofar <- function(x,rank=NULL){
  fit$Y - fitted.sofar(x=x,rank=rank)
}
summary.sofar <- function(fit,rank=NULL){
  yhat <- fitted.sofar.Z(fit,rank)
  totalVar <- function(x){
    sum(diag(crossprod(x)))
  }
  totalV <- totalVar(fit$Y)
  R2 <- 1-totalVar(fit$Y-yhat)/totalV
  cat("rank=",fit$rank,"\n")
  cat("R2=",R2,"\n")
  cp <- cumsum(fit$D)
  names(cp) <- paste0("factor1-",1:length(cp))
  cp
}

dev_ <- function(){
  rm(list=ls())
  load(file="Lee_SOFAR.rData")
  source("sofarUtils.R")
  table(group.sofar(fit,rank=3))
  f <- fitted(fit,rank=3)
  e <- residuals(fit,rank=3)
  summary(fit,rank=3)
}