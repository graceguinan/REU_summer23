
# full code in bclust handout

rm(list=ls())
library("dplyr")


# load data
## common signature quant SOFAR run
X <- read.csv("data/HRMS_data.csv")

# select non-zero columns (common signature)
X <- X %>% select_if(~min(., na.rm = TRUE) > 0)
X1 <- as.matrix(X[,-c(1:4)])
dim(X1)


# normalize data
## log transformation
X2 <- log(X1,base=10)

## center and scale
alpha = .05
center <- apply(X2,2,mean)
X2b <- sweep(X2, 2, center, FUN = "-")
scale <- apply(X2b,2, function(x) diff(quantile(x,probs=c(alpha,1-alpha))))
X3 <- sweep(X2b, 2, scale, FUN = "/")


# run SOFAR
if(!interactive()){
  library(rrpack)
  X3 = as.matrix(X3)
  fit <- sofar(Y = X3 ,X=diag(nrow(X3)),nrank = 12,ic.type="GIC",control =
                 list(nlam=100,lam.min.factor=1e-7))
  save(fit, center, scale, file = "fit_common_signature.rData") 
}

## load fit + files
library(RColorBrewer)
library(ggplot2)

load(file="fit/fit_common_signature.rData")
source("helper_functions/sofarUtils.R")
source("helper_functions/image_sofar_3.R")
summary(fit)

## percent explained by each layer
r2 = 86
svd1 = fit[["D"]]
layers = svd1^2/sum(svd1^2)
layers * r2


# plot
## fitted values
yhat <- fitted.sofar.Z.rank(fit,rank=1)
image.sofar(yhat, legend.title = "")+guides(fill = guide_legend(reverse=TRUE))

## y cutoff
min(which(yhat[,1]<0))

## x cutoff
max(which(yhat[1,]>0))
min(which(yhat[1,]<0))

## plot with annotations
yhat <- fitted.sofar.Z.rank(fit,rank=1)
image.sofar(yhat, legend.title = "")+guides(fill = guide_legend(reverse=TRUE))+
  annotate("segment", x = 810, xend = 810, y = 0, yend = 145,colour = "black", linewidth = 1, alpha = .7)+ 
  annotate("segment", x = 943, xend = 943, y = 0, yend = 145,colour = "black", linewidth = 1, alpha = .7)+ 
  annotate("segment", x = 0, xend = 2201, y = 68, yend = 68, colour = "black", linewidth = 1, alpha = .7)+
  annotate("text", x = -100, y = 104.5, label = "Aloha Samples", angle = 90, size = 4)+
  annotate("text", x = -100, y = 34, label = "Bats Samples", angle = 90, size = 4)+
  annotate("text", x =410 , y = -3, label = "Group 1", size = 4)+
  annotate("text", x =876 , y = -3, label = "Group 2", size = 4)+
  annotate("text", x =1575 , y = -3, label = "Group 3", size = 4)+
  ggtitle("Layer 1")+ylab("")+xlab("")+
  coord_cartesian(xlim = c(0, 2201), ylim = c(0,145),  clip = 'off')+
  theme(axis.text.x = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_blank(), axis.title=element_text(size=45))