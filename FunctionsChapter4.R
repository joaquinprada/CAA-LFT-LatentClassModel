## Generate a truncated distribution - from Keith Goldfeld in R-bloggers
rnormt <- function(n, range, mu, s = 1) {
  
  # range is a vector of two values
  
  F.a <- pnorm(min(range), mean = mu, sd = s)
  F.b <- pnorm(max(range), mean = mu, sd = s)
  
  u <- runif(n, min = F.a, max = F.b)
  
  qnorm(u, mean = mu, sd = s)
}

modelSensSpecCAA <- function(nInd,ID){
  nInf <- rbinom(1,nInd,prob = P[ID])
  Infectlvl <- rgamma(nInf,shape = sh[ID], rate = rt[ID])
  
  Intensity <- c(rnormt((nInd-nInf),range = c(0,10), mu = 0, s = sdCAA[ID]),
                 rnorm(nInf, mean = 10 / (1 + exp(-param1[ID]*(Infectlvl-param2[ID]))), sd = sdCAA[ID]))
  
  sens <- array(NA,11)
  spec <- array(NA,11)
  for (t in 0:10){
    sens[t+1] <- length(which(Intensity[(nInd-nInf+1):nInd]>=t))/nInf
    spec[t+1] <- length(which(Intensity[1:(nInd-nInf)]<t))/(nInd-nInf)
  }
  return(cbind(sens,spec))
}

modelSensSpecCCA <- function(nInd,ID){
  nInf <- rbinom(1,nInd,prob = P[ID])
  Infectlvl <- rgamma(nInf,shape = sh[ID], rate = rt[ID])
  
  Intensity <- c(rnormt((nInd-nInf),range = c(0,9), mu = 0, s = sdCCA[ID]),
                 rnorm(nInf, mean = 9 / (1 + exp(-multim[1]*(Infectlvl-multim[2]))), sd = sdCCA[ID]))
  
  sens <- array(NA,10)
  spec <- array(NA,10)
  for (t in 0:9){
    sens[t+1] <- length(which(Intensity[(nInd-nInf+1):nInd]>=t))/nInf
    spec[t+1] <- length(which(Intensity[1:(nInd-nInf)]<t))/(nInd-nInf)
  }
  return(cbind(sens,spec))
}


