####################
###    dWOLS     ###
### single stage ###
####################

library(parallel)
library(foreach)
library(DTRreg)
library(doParallel)

rm(list=ls()) # avoid confusion due to same variable name
if (getwd() == '/Users/shamalama'){
  setwd('~/Documents/McGill/Research/Thesis/simulations/')
}

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

psi <- c(-0.5, 1)  # (psi0, psi1)

num.iter <- c(10, 100)
num.patients <- 1000

use.dtrreg <- TRUE # if analysis will be performed alongside DTRreg package

create.data <- function(n){
  X <- rnorm(n=n)
  A <- rbinom(n, 1, 1/(1+exp(-X)))
  Y <- rnorm(n=n, mean=exp(X) - X^3 + A*(psi[1] + X*psi[2]), sd=1)
  
  return (cbind(Y, A, X))
}


psi.func <- function(.A, .Ahat, .H.beta, .H.psi, .Y){
  n <- length(.A)
  
  W   <- diag(abs(.A - .Ahat))
  inv <- solve((t(.H.beta)) %*% W %*% .H.beta)
  
  Delta <- diag(n) - .H.beta %*% inv %*% t(.H.beta) %*% W
  Cross <- t(.H.psi) %*% W %*% diag(.A) %*% Delta
  
  psi <- solve(Cross %*% diag(.A) %*% .H.psi) %*% Cross %*% .Y
  
  return (as.vector(psi))
}

psi.func.lm <- function(.A, .Ahat, .H.beta, .H.psi, .Y){
  n <- length(.A)
  
  model <- lm(.Y ~ 0+cbind(.H.beta, .A*.H.psi), weights=abs(.A - .Ahat))
  coefs <- as.vector(model$coefficients)
  
  stopifnot(length(coefs) == (dim(.H.beta)[2] + dim(.H.psi)[2]))
  
  output <- coefs[dim(.H.beta)[2]+1:length(coefs)]
  
  return (output)
}

for (n in num.iter){

  
  f <- function(i){
    temp.data <- create.data(num.patients)
    Y <- temp.data[,1]
    A <- temp.data[,2]
    X <- temp.data[,3]
    
    Ahat.inc <- as.vector(fitted(glm(A ~ 1, family=binomial)))
    Ahat.cor <- as.vector(fitted(glm(A ~ X, family=binomial)))
    
    x.beta.inc <- cbind(1, X)
    x.beta.cor <- cbind(1, exp(X), X^3)
    x.psi <- cbind(1, X)
    
    
    # both models incorrect
    psi.model1   <- psi.func.lm(A, Ahat.inc, x.beta.inc, x.psi, Y)
    
    # treatment model incorrect
    psi.model2   <- psi.func.lm(A, Ahat.inc, x.beta.cor, x.psi, Y)
    
    # outcome model incorrect
    psi.model3   <- psi.func.lm(A, Ahat.cor, x.beta.inc, x.psi, Y)
    
    # both models correct
    psi.model4   <- psi.func.lm(A, Ahat.cor, x.beta.cor, x.psi, Y)
    
    return (rbind(psi.model1,
                  psi.model2,
                  psi.model3,
                  psi.model4))
  }
  ptm <- proc.time() # start timer
  obj <- foreach(x = 1:n, .combine=rbind) %dopar% f(x)
  
  indices <- seq(1, n)*4
  dwols.model1.psi <- obj[indices-3,]
  dwols.model2.psi <- obj[indices-2,]
  dwols.model3.psi <- obj[indices-1,]
  dwols.model4.psi <- obj[indices,]
  
  
  png(paste('dwols-single-stage-psi_n', n, '.png', sep=''), width=10, height=6, units='in', res=300)
  par(mfrow=c(1, 2))
  boxplot(cbind(dwols.model1.psi[,1],
                dwols.model2.psi[,1],
                dwols.model3.psi[,1],
                dwols.model4.psi[,1]),
          ylim=c(-1, 0.8))
  abline(h=psi[1], lty=2, lwd=1)
  title(xlab='Model type',
        ylab=expression(paste(psi[0], '')),
        main=expression(paste(psi[0], ' values')))
  
  boxplot(cbind(dwols.model1.psi[,2],
                dwols.model2.psi[,2],
                dwols.model3.psi[,2],
                dwols.model4.psi[,2]),
          ylim=c(-0.8, 4))
  abline(h=psi[2], lty=2, lwd=1)
  title(xlab='Model type',
        ylab=expression(paste(psi[1], '')),
        main=expression(paste(psi[1], ' values')))
  dev.off()
}
stopCluster(cl)
# stop timer
proc.time() - ptm
