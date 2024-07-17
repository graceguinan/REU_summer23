
# full code in sofar handout

# packages and data
rm(list=ls())
library(rrpack)
library(stringr)
library(ggplot2)
library(dplyr)

## load data
X <- read.csv("data/HRMS_data.csv")
Y <- read.csv("data/FDOM_data.csv")


# normalization
## center and scale (quant)
X2<-sapply(X[,-c(1:3)], function(x) quantile(x,probs = c(.05,.95) ))
X3 <-as.data.frame(t(X2))

## make new column subtracting 95th and 5th percentile
X3$Sj <- X3$`95%` - X3$`5%`

## select only nonzero values
X4 <- X3[X3$Sj != 0,]

## delete these compounds from original dataset
X4$ID <- row.names(X4)
XX <- X[,c("Sample", "Location", "Depth",X4$ID)]

## add Sj into X dataset as a row
Sv<-X4$Sj
Sv<-t(c(NA, NA, NA, Sv))
colnames(Sv)<- colnames(XX)
XX<- rbind(Sv, XX)

## make Z dataset (divide X by Sj)
XX2<- XX[-c(1:3)]
Z <- XX2/XX2[rep(1,nrow(XX2)),]
Z<- cbind(XX[,1:3],Z)
Z<- Z[-1,]


# merge x and y
idlist <- list()
for (i in 1:145){
  q<-str_sub(Z$Sample[i], -1, -1)
  idlist[[length(idlist)+1]] <-q
}
Z$Profile <-idlist

## merge data together
Z2 <- merge(Z,Y[,c("Profile","Depth","Location")],by=c("Profile","Depth","Location"))
y3 <- merge(Y,Z2[,c("Profile","Depth","Location")],by=c("Profile","Depth","Location"))


# SOFAR
if(!interactive()){
  z31 <- as.matrix(Z2[,-c(1:4)])
  y31 <- as.matrix(y3[,-c(1:3)])
  fit_quantile_all1234<- sofar(y31, z31, nrank = 3,
                               control=list(nlam = 100, 
                                            lam.min.factor=1e-10,
                                            lam.max.factor=5),ic.type="BIC")
  
  save(fit_quantile_all1234,file="fit/sofarFit2ALL_quantile1234.rData")}

## load fit
load(file="fit/sofarFit2ALL_quantile1234.rData")



# fitted values
fitted_sofar <- function(x){
  with(x,X%*%U%*%diag(D)%*%t(V))
}
yhat <- fitted_sofar(fit_quantile_all1234)

## R-square
totalVar <- function(x){
  sum(diag(crossprod(x)))
}
1-totalVar(fit_quantile_all1234$Y-yhat)/totalVar(fit_quantile_all1234$Y)

## plot fitted versus observed
linkY <- y3[c(1:3)]
linkY$i <- as.numeric(row.names(linkY))

data_ <- data.frame(
  y=c(fit_quantile_all1234$Y),
  fit=c(yhat),
  i=c(row(fit_quantile_all1234$Y)),
  j=c(col(fit_quantile_all1234$Y))
)
data_1 <- full_join(linkY, data_, by = "i")
avgdat <- aggregate(cbind(y, fit)~j +Location +Depth, FUN= mean, data = data_1)

## graph with location and depth
ggplot(data =avgdat,aes(y = y, x = Depth, group = Location))+geom_line(linetype=1, aes(color = Location))+
  geom_line(aes(x=Depth,y=fit, color=Location),linetype=2)+
  facet_wrap(~j,scales = "free")+coord_flip() + scale_x_reverse()+ylim(0,1.3)+scale_color_manual(values = c("red", "blue"))+theme(axis.title.x =element_blank())

## scatter plot
ggplot(avgdat,aes(x=y,y=fit, color = Location))+geom_point()+
  facet_wrap(~j,scales="free")+xlab("observed")+ylab("fitted")+scale_color_manual(values = c("red", "blue"))


# Venn Diagram
#install.packages("VennDiagram")
#install.packages("gplots")
library(VennDiagram)
library(gplots)

## make T/F dataset for if compounds are inclued in each grouping
ULL2 <- as.data.frame(fit_quantile_all1234[["U"]])
ULL2[ULL2 >= .01] <-1
ULL2[ULL2 <= -.01] <-1
ULL2$V1 <- ULL2$V1==1
ULL2$V2 <- ULL2$V2==1
ULL2$V3 <- ULL2$V3==1

## Create logical vectors
set1 <- ULL2$V1
set2 <- ULL2$V2
set3 <- ULL2$V3

venn_data <- list(
  col1 = which(set1),
  col2 = which(set2),
  col3 = which(set3)
)
unselected_count <- sum(!set1 & !set2 & !set3)

## Generate the Venn diagram
venn.plot <- venn.diagram(
  x = venn_data,
  category.names = c("1,2,4", "3", "1,2,4"),
  output = FALSE,
  main = "Count of componds in each group", 
  filename = NULL,
  fill = c("red", "darkorange", "blue"),   
  alpha = 0.5,                        
  cat.col = c("red", "darkorange", "blue"), 
  cat.cex = 1.5,                     
  cat.fontface = "bold"            
)
grid.newpage()
grid::grid.draw(venn.plot)
grid.text(paste(unselected_count), x = 0.9, y = 0.1, gp = gpar(fontsize = 12, col = "black"))


# heatmap
library(viridis)
library(cowplot)

## load metadata and run PCA
metadat <- read.csv("data/metadat.csv")
meta.pca <-princomp(metadat[-1], cor = TRUE)
ggdata_ <- data.frame(meta.pca$scores, ID = metadat$ID)

## make heat map
source("helper_functions/heatmap.R")
ID_quant <- as.list(read.csv("data/quant_ids.csv")$ID)
heatmap(fit_quantile_all1234, ID_quant, ggdata_)

## make PCA
loads <- as.data.frame(meta.pca$loadings[,c(1:2)])
pca_heatmap(ggdata_, loads)