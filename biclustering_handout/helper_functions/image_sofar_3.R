image.sofar.mesh <- function(mat){
  ## matrix to data frame.
  data.frame(y=c(row(mat)),x=c(col(mat)),z=c(mat))
}
image.sofar <- function(x,brk=NULL,legend.label=NULL,
                           pal="YlOrRd",zlim=NULL,
                           legend.title="unit",xlab="feature",ylab="sample"){
  if(is.matrix(x)){
    mesh <- image.sofar.mesh(x)
  }
  else if(class(x)=="list"){
    meshfit <- image.sofar.mesh(x[[1]])
    meshobs <- image.sofar.mesh(x[[2]])
    mesh <- rbind(cbind(meshfit,type="fit"),cbind(meshobs,type="data"))
  }
  else{
    stop("x must be matrix of list.\n")
  }

  ## legend break
  if(is.null(zlim)){
    if(is.null(brk)){
      brk <- pretty(mesh$z)
    }
  }else{
    if(is.null(brk)){
      brk <- pretty(c(zlim,mesh$z))
    }else{
      stopifnot(max(brk)>max(zlim))
    }
  }
  
  mesh$col <- cut(mesh$z,breaks=brk,labels=legend.label)
  
  col_ <- brewer.pal(n=length(brk)-1,name = pal)
  #browser()
  base <- ggplot(data=mesh,aes(x=x,y=y,fill=col))+
    geom_raster(interpolate = FALSE)+xlab(xlab)+ylab(ylab)+theme_bw()+
    scale_fill_manual(name=legend.title,values=col_)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))
  
  if(is.list(x)){
    return(base+facet_grid(type~.))
  }else{
    return(base)
  }
}

dev_ <- function(){
  rm(list=ls())
  library(ggplot2)
  library(RColorBrewer)
  source("sofarUtils.R")
  source("image_sofar_3.R")
  load(file="Lee_SOFAR_C1Q1M5.rData")
  ## load masses with common signature
  load(file="HRMS_Clst.rData")
  min_ <- apply(X,2,min)
  X1 <- X[,min_>0]
  image.sofar(X1)
  quantile(c(X1/1e6))
  image.sofar(X1/1e6,brk=c(0,10,30,100,1000,10000),zlim=9999,
              legend.label = c("<10","11-30","30-100",
                               "100-1000","1000-6700"))
  summary.sofar(fit)
  g <- group.sofar(fit,rank=4)

  image.sofar(X1[,order(g)]/1e6,brk=c(0,10,30,100,1000,10000),zlim=9999,
              legend.label = c("<10","11-30","30-100",
                               "100-1000","1000-6700"))
  
  debug(fitted.sofar)
  fit <- fitted(fit,rank=4,log=T,center = center,scale=scale_,base = 10)
  
  image.sofar(fit[,order(g)]/1e6,brk=c(0,10,30,100,1000,10000),zlim=9999,
              legend.label = c("<10","11-30","30-100",
                               "100-1000","1000-6700"))
  
  x <- list(fit[,order(g)]/1e6,X1[,order(g)]/1e6)
  image.sofar(x,brk=c(0,10,30,100,1000,10000),zlim=9999,
              legend.label = c("<10","11-30","30-100",
                               "100-1000","1000-6700"))
  
  ## Figure 1 of Lee et al. 2010
  #image.sofar(fit,col=brewer.pal(n=9,name="RdBu"),ptype = "data")

}