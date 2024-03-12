library(rstudioapi)
library(runjags)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory

## Load data
dt <- read.csv("ch4_bug22_including_CAA_concentration.csv")

## Prepare data
KK <- as.matrix(dt[,5:10])
CCA <- dt$POC.CCA-1
CAA <- dt$CAA_GNP3
CAA[which(CAA==1)]<-0
CAA[which(CAA==2)]<-1
CAA[which(CAA==3)]<-2
CAA[which(CAA==3.5)]<-3

## Load posterior CCA logistic function
dtprior <- read.csv("kints.csv")

## Fit a multivariate normal distribution for the two parameters
library(MGMM)
fitParams <- FitGMM(data = as.matrix(dtprior[,c(4,5)]))

multim <- fitParams@Mean
multicov <- fitParams@Covariance


## Set seed ##
.RNG.seed <- function(chain)
  return( switch(chain, "1"= 1, "2"= 2) )
.RNG.name <- function(chain)
  return( switch(chain, "1" = "base::Super-Duper", "2" = "base::Wichmann-Hill") )

## Inits ##
N = nrow(dt)
nKK = ncol(KK)
status = rep(1,N)
param1 = 1
param2 = 0
CAAIntensity <- array(rep(CAA,2),dim = c(N,2))
CCAIntensity <- array(rep(CCA,2),dim = c(N,2))

m <- "model{

  for (n in 1:N){
    ## Status
    status[n] ~ dbern(P)
    
    infectIntensity[n] ~ dgamma(sh,rt)
    
    ## KK 
    for (r in 1:nKK){
      KK[n,r] ~ dnegbin(k/(infectIntensity[n]*status[n]+k),k)
    }
    
    ## CCA component
    CCAIntensity[n,1] ~ dnorm(0,tau)T(0,9)
    CCAIntensity[n,2] ~ dnorm(9 / (1 + exp(-multiParam[1]*(infectIntensity[n]-multiParam[2]))),tau)
    
    CCA[n] ~ dround(CCAIntensity[n,status[n]+1],0)
    
    ## CAA component
    CAAIntensity[n,1] ~ dnorm(0,tauCAA)
    CAAIntensity[n,2] ~ dnorm(10 / (1 + exp(-param1*(infectIntensity[n]-param2))),tauCAA)
    
    CAA[n] ~ dround(CAAIntensity[n,status[n]+1],0)
    ## CAA component
    #CAA[n] ~ dnorm(10 / (1 + exp(-param1*(infectIntensity[n]*status[n]-param2))),tauCAA)T(0,)
    
  }
  
  ## Prior 
  P ~ dbeta(1,1)
  sh ~ dgamma(83.88592,125.70331)
  rt1 ~ dbeta(47.13542,633.08366)
  rt <- rt1/(1-rt1)
  k ~ dgamma(0.001,0.001)
  tau ~ dgamma(0.001,0.001)
  tauCAA ~ dgamma(0.001,0.001)
  multiParam[1:2] ~ dmnorm.vcov(multim, multicov)
  param1 ~ dgamma(0.001,0.001)
  param2 ~ dnorm(0,10^-6)
  
  #inits# .RNG.seed, .RNG.name, status, param1, param2, CCAIntensity, CAAIntensity
  #data# N, nKK, KK, CCA, CAA, multim, multicov
  #monitor# P, sh, rt, k, tau, tauCAA, param1, param2

}"

Results <- run.jags(m, burnin=1000, sample=5000, thin=1, n.chains=2, jags.refresh = 1, method = 'parallel',
                    plots = F, silent.jags = F)

plot(Results)

#save.image("Results.RData")

###################################################################
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path)) #set file location as working directory
load("Results.RData")
source("FunctionsChapter4.R")

## Set draws
set.seed(100)
ndraws <- 1000
IDs <- sample.int(10000,ndraws)

## Calc Sens/Spec
P <- c(Results$mcmc[[1]][,"P"],Results$mcmc[[2]][,"P"])[IDs]
sh <- c(Results$mcmc[[1]][,"sh"],Results$mcmc[[2]][,"sh"])[IDs]
rt <- c(Results$mcmc[[1]][,"rt"],Results$mcmc[[2]][,"rt"])[IDs]
sdCAA <- sqrt(1/c(Results$mcmc[[1]][,"tauCAA"],Results$mcmc[[2]][,"tauCAA"]))[IDs]
sdCCA <- sqrt(1/c(Results$mcmc[[1]][,"tau"],Results$mcmc[[2]][,"tau"]))[IDs]
param1 <- c(Results$mcmc[[1]][,"param1"],Results$mcmc[[2]][,"param1"])[IDs]
param2 <- c(Results$mcmc[[1]][,"param2"],Results$mcmc[[2]][,"param2"])[IDs]

sensSpecEstimatesCAA <- sapply(1:length(IDs),function(x){modelSensSpecCAA(1000,x)})
QsenspcCAA <- apply(sensSpecEstimatesCAA,1,quantile,c(.025,.5,.975))

sensSpecEstimatesCCA <- sapply(1:length(IDs),function(x){modelSensSpecCCA(1000,x)})
QsenspcCCA <- apply(sensSpecEstimatesCCA,1,quantile,c(.025,.5,.975))


pdf("ROCcurve.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(NA,xlim = c(0,1),ylim = c(0,1),axes = F, xlab = "", ylab = "")
## CAA
points(1-QsenspcCAA[2,12:22],QsenspcCAA[2,1:11], pch=19)
lines(1-QsenspcCAA[2,12:22],QsenspcCAA[2,1:11])
arrows(x0=1-QsenspcCAA[2,12:22],y0=QsenspcCAA[1,1:11],
       y1=QsenspcCAA[3,1:11], angle=90, code = 3,length = .05)
arrows(x0=1-QsenspcCAA[1,12:22],y0=QsenspcCAA[2,1:11],
       x1=1-QsenspcCAA[3,12:22], angle=90, code = 3,length = .05)

## CCA
points(1-QsenspcCCA[2,11:20],QsenspcCCA[2,1:10], pch=19, col="darkgreen")
lines(1-QsenspcCCA[2,11:20],QsenspcCCA[2,1:10], col="darkgreen")
arrows(x0=1-QsenspcCCA[2,11:20],y0=QsenspcCCA[1,1:10],
       y1=QsenspcCCA[3,1:10], angle=90, code = 3,length = .05, col="darkgreen")
arrows(x0=1-QsenspcCCA[1,11:20],y0=QsenspcCCA[2,1:10],
       x1=1-QsenspcCCA[3,11:20], angle=90, code = 3,length = .05, col="darkgreen")

axis(1)
axis(2)

mtext("Sensitivity",side=2,cex=1,line=1.2)
mtext("1 - Specificify",side=1,cex=1,line=1.2)

legend("bottomright",c("CAA","CCA"),
       col=c("black","darkgreen"),
       bty='n',cex=0.75,lty=1)

dev.off()


###### AUC with 50% quantile
library(sf)
###CAA
polygon <- st_sfc(st_polygon(list(cbind(c(1-QsenspcCAA[2,12:22],0,1,1-QsenspcCAA[2,12]),
                                        c(QsenspcCAA[2,1:11],0,0,QsenspcCAA[2,1])))))

plot(polygon,axes = TRUE)
st_area(polygon)

###CCA
polygon <- st_sfc(st_polygon(list(cbind(c(1-QsenspcCCA[2,11:20],0,1,1-QsenspcCCA[2,11]),
                                        c(QsenspcCCA[2,1:10],0,0,QsenspcCCA[2,1])))))

plot(polygon,axes = TRUE)
st_area(polygon)


###### CAA Curve
x <- 0:100
library(scales)

meanC <- rowMeans(sapply(1:length(IDs),function(i){10 / (1 + exp(-param1[i]*(x-param2[i])))}))

pdf("CAAcurve.pdf",width=7,height=4.64)
par(font=2, cex.axis=0.75, lwd=2, mar=c(2.1,2.2,1.2,0)+0.1,mgp=c(3,0.4,0))

plot(NA,xlim = c(0,100),ylim = c(0,10),axes = F, xlab = "", ylab = "")

for (i in 1:length(IDs)){
  lines(x,10 / (1 + exp(-param1[i]*(x-param2[i]))), col=alpha("Grey",.2))
}

lines(x,meanC)

axis(1, at = seq(0,100,by=20), labels = seq(0,100,by=20)*24)
axis(2, at = 0:10, labels = c(1:3,3.5,4:10))

mtext("CAA LFT G-score",side=2,cex=1,line=1.2)
mtext("Intensity of Infection (EPG)",side=1,cex=1,line=1.2)

dev.off()