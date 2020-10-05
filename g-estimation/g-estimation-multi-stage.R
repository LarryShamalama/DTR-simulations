######################
###  G-estimation  ###
###  Three stage   ###
######################

# takes around 7 minutes to run

library(parallel)
library(foreach)
library(DTRreg)
library(doParallel)

rm(list=ls())
if (getwd() == '/Users/shamalama'){
  setwd('~/Documents/McGill/Research/Thesis/simulations/multi-stage/G-estimation')
}
ptm <- proc.time() # start timer

use.dtrreg <- FALSE # if analysis will be performed alongside DTRreg package

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)

# Vary setups
# need to make sure that psi0j + psi1j*Xj > 0
# is a rule that makes sense for all stages j

# all same for simplicity reasons
psi1 <- c(-0.5, 1)
psi2 <- c(-0.5, 1)
psi3 <- c(-0.5, 1)


num.iter     <- c(100, 1000, 10000)
num.stages   <- 3
n <- 1000 # number of patients
  
psi.func <- function(.A, .Ahat, .H.beta, .H.psi, .Y){
  W   <- diag(.A - .Ahat) - (.A - .Ahat)*.H.beta %*% solve(t(.H.beta) %*% .H.beta) %*% t(.H.beta)
  psi <- as.vector(solve( t(.H.psi) %*% W %*% (.A*.H.psi)) %*% t(.H.psi) %*% W %*% .Y)  
  
  return (psi)
}

for (iter in num.iter){
  
  
  f <- function(i){
    psi1 <- c(-0.5, 1)
    psi2 <- c(-0.5, 1)
    psi3 <- c(-0.5, 1)
    
    # set-up (data generation)
    
    X1 <- rnorm(n)
    A1 <- rbinom(n, 1, 1/(1+exp(-X1)))
    
    X2 <- A1 + rnorm(n)
    A2 <- rbinom(n, 1, 1/(1+exp(-X2)))
    
    X3 <- A2 + rnorm(n)
    A3 <- rbinom(n, 1, 1/(1+exp(-X3)))
    
    A1.opt <- as.numeric(psi1[1]+psi1[2]*X1>0)
    A2.opt <- as.numeric(psi2[1]+psi2[2]*X2>0)
    A3.opt <- as.numeric(psi3[1]+psi3[2]*X3>0)
        
    Y  <- 1 + X1 + A1*(psi1[1]+psi1[2]*X1) - abs(A2-A2.opt)*abs(psi2[1]+psi2[2]*X2) - abs(A3-A3.opt)*abs(psi3[1]+psi3[2]*X3)
    
    
    # stage 3 estimation
    
    A3hat.cor <- as.vector(fitted(glm(A3 ~ X3, binomial)))
  
    Y3 <- Y
    H3.beta.inc <- cbind(rep(1, n))
    H3.psi      <- cbind(rep(1, n), X3)
  
    # treatment-free model incorrect
    psi.hat.3 <- psi.func(A3, A3hat.cor, H3.beta.inc, H3.psi, Y3)
    
    A3opt.hat <- as.integer(H3.psi %*% psi.hat.3 > 0)

    
    # stage 2 estimation
    
    A2hat.cor <- as.vector(fitted(glm(A2 ~ X2, binomial)))
    
    Y2 <- Y3 + (A3 != A3opt.hat)*(abs(H3.psi %*% psi.hat.3))
    H2.beta.inc <- cbind(rep(1, n))
    H2.psi      <- cbind(rep(1, n), X2)
    
    # treatment-free model incorrect
    psi.hat.2 <- psi.func(A2, A2hat.cor, H2.beta.inc, H2.psi, Y2)
    
    A2opt.hat <- as.integer(H2.psi %*% psi.hat.2 > 0)
    
    
    # stage 1 estimation
    
    A1hat.inc <- as.vector(fitted(glm(A1 ~ 1, binomial)))
    
    Y1 <- Y2 + (A2 != A2opt.hat)*(abs(H2.psi %*% psi.hat.2))
    H1.beta.cor <- cbind(rep(1, n), X1)
    H1.psi      <- cbind(rep(1, n), X1)
    
    psi.hat.1 <- psi.func(A1, A1hat.inc, H1.beta.cor, H1.psi, Y1)
    
    A1opt.hat <- as.integer(H2.psi %*% psi.hat.1 > 0)
    
    blip.mod  <- list(~X1, ~X2, ~X3) # has to be correct
    treat.mod <- list(A1~X1, A2~X2, A3~X3)
    tf.mod    <- list(~X1, ~1, ~1)
    
    data <- data.frame(X1, A1, X2, A2, X3, A3)
    
    # Using DTRreg
    
    
    return (rbind(psi.hat.1, psi.hat.2, psi.hat.3))
  }
  
  ptm <- proc.time() # start timer
  obj <- foreach(x = 1:iter, .combine=rbind) %dopar% f(x)
  
  indices <- seq(1, iter)*3
  psi.acc.1 <- obj[indices-2,]
  psi.acc.2 <- obj[indices-1,]
  psi.acc.3 <- obj[indices,]
  
  
  ################
  ### BOXPLOTS ###
  ################
  
  # Stage 3 psi0
  
  png(paste('multi-stage-psi_n', iter, '.png', sep=''), width=10, height=6, units='in', res=300)
  par(mfrow=c(1, 2))
  boxplot(cbind(psi.acc.3[,1],
                psi.acc.2[,1],
                psi.acc.1[,1]),
          names=c(expression(psi[30]),
                  expression(psi[20]),
                  expression(psi[10])))
  abline(h=psi3[1], lty=2)
  title(main=expression(paste(psi[0], ' values')))
  
  boxplot(cbind(psi.acc.3[,2],
                psi.acc.2[,2],
                psi.acc.1[,2]),
          names=c(expression(psi[31]),
                  expression(psi[21]),
                  expression(psi[11])))
  abline(h=psi3[2], lty=2)
  title(main=expression(paste(psi[1], ' values')))
  dev.off()
}

# stop timer
proc.time() - ptm
