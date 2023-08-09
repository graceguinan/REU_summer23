# SOFAR on ALL data 
#LOG transformed

library(readr)
library(stringr)
library(dplyr)
library(rrpack)
# set working directory

#load datasets, name X and Y
X <- read_csv("HRMS and FDOM data BATS ALOHA.xlsx - HRMS data.csv")
Y <- read_csv("HRMS and FDOM data BATS ALOHA.xlsx - FDOM data.csv")

# delete component 5 (fmax5) from Y data set
Y <- Y[,-8]


# make id number column in HRMS data
X2 <- X
v <- as.character(X2[1,])
names(X2) <- v
X2 <- X2[-1,]
colnames(X2)[1] <- "Sample"
colnames(X2)[2] <- "Location"
colnames(X2)[3] <- "Depth"

# correct names of bats 3800 and aloha 400
X2[57,1] <- "BATS_3800m_3"
X2[80,1] <- "HOT_400m_3"

idlist <- list()
for (i in 1:145){
  q<-str_sub(X2$Sample[i], -1, -1)
  idlist[[length(idlist)+1]] <-q
}
X2$Profile <-idlist
X2 <-X2[,c(6153, 1:6152)]
colnames(Y)[3] <- "Depth"

# merge data together
x3 <- merge(X2,Y[,c("Profile","Depth","Location")],by=c("Profile","Depth","Location"))
y3 <- merge(Y,x3[,c("Profile","Depth","Location")],by=c("Profile","Depth","Location"))

# edit and log transform
x31 <- x3[,-c(1:4)]
y31 <- y3[,-c(1:3)]
x311 <- log(x31, base= 10)+1
x311[x31 == 0] <-0

#run SOFAR
y31 = as.matrix(y31)
x311 = as.matrix(x311)
fit_log_ALL_1234 <- sofar(y31, x311, nrank = 3,control=list(nlam = 100, lam.min.factor=1e-10,lam.max.factor=5),ic.type="BIC")

save(fit_log_ALL_1234,file="sofarFitALL_log1234.rData")
