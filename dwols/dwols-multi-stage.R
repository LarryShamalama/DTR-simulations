######################
###     dWOLS      ###
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

use.dtrreg <- TRUE # if analysis will be performed alongside DTRreg package

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


num.iter     <- c(10)
num.stages   <- 3
n <- 1000 # number of patients

psi.func <- function(.A, .Ahat, .H.beta, .H.psi, .Y){
  W   <- diag(abs(.A - .Ahat))
  inv <- solve((t(.H.beta)) %*% W %*% .H.beta)
  one <- (t(.A*.H.psi)) %*% W %*% .Y
  B   <- (t(.A*.H.psi)) %*% W %*% (.A*.H.psi)
  two <- (t(.A*.H.psi)) %*% W %*% .H.beta %*% inv %*% t(.H.beta) %*% W %*% .Y
  D   <- (t(.A*.H.psi)) %*% W %*% .H.beta %*% inv %*% t(.H.beta) %*% W %*% (.A*.H.psi)
  
  psi <- solve(B - D) %*% (one - two)
  
  return (as.vector(psi))
}

for (iter in num.iter){
  
  
  f <- function(i, use.package){
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
    
    invX3 <- 1/X3
    Y3 <- Y
    H3.beta.inc <- cbind(rep(1, n), invX3)
    H3.psi      <- cbind(rep(1, n), X3)
    
    # treatment- free model incorrect
    psi.hat.3 <- psi.func(A3, A3hat.cor, H3.beta.inc, H3.psi, Y3)
    
    A3opt.hat <- as.integer(H3.psi %*% psi.hat.3 > 0)
    
    
    # stage 2 estimation
    
    A2hat.cor <- as.vector(fitted(glm(A2 ~ X2, binomial)))
    
    invX2 <- 1/X2
    Y2 <- Y3 + (A3 != A3opt.hat)*(abs(H3.psi %*% psi.hat.3))
    H2.beta.inc <- cbind(rep(1, n), invX2)
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
    
    if (!use.package){
      return (rbind(psi.hat.1, psi.hat.2, psi.hat.3))
    }
    else{

      blip.mod  <- list(~X1, ~X2, ~X3) # has to be correct
      treat.mod <- list(A1~1, A2~X2, A3~X3)
      tf.mod    <- list(~X1, ~invX2, ~invX3)
      
      data <- data.frame(X1, A1, X2, A2, X3, A3, Y)
      
      model <- DTRreg(Y, blip.mod, treat.mod, tf.mod, method='dwols')
      
      return (rbind(psi.hat.1, 
                    psi.hat.2, 
                    psi.hat.3,
                    model$psi[[1]], 
                    model$psi[[2]], 
                    model$psi[[3]]))
    }
  }
  
  
  if (!use.dtrreg){
    ptm <- proc.time() # start timer
    obj <- foreach(x = 1:iter, .combine=rbind) %dopar% function(x){f(x, use.dtrreg)}
    
    indices <- seq(1, iter)*3
    psi.acc.1 <- obj[indices-2,]
    psi.acc.2 <- obj[indices-1,]
    psi.acc.3 <- obj[indices,]
    
    
    ################
    ### BOXPLOTS ###
    ################
    
    png(paste('dwols-multi-stage-psi_n', iter, '.png', sep=''), width=10, height=6, units='in', res=300)
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
  
  else{
    ptm <- proc.time() # start timer
    obj <- foreach(x = 1:iter, .combine=rbind) %dopar% function(x){return (f(x, use.dtrreg))}
    
    indices <- seq(1, iter)*6
    psi.acc.1 <- obj[indices-5,]
    psi.acc.2 <- obj[indices-4,]
    psi.acc.3 <- obj[indices-3,]
    psi.acc.1.dtrreg <- obj[indices-2,]
    psi.acc.2.dtrreg <- obj[indices-1,]
    psi.acc.3.dtrreg <- obj[indices,]
    
    png(paste('dwols-multi-stage-dtrreg-psi0_n', iter, '.png', sep=''), width=10, height=6, units='in', res=300)
    par(mfrow=c(1, 2))
    boxplot(cbind(psi.acc.3[,1],
                  psi.acc.2[,1],
                  psi.acc.1[,1]),
            names=c(expression(psi[30]),
                    expression(psi[20]),
                    expression(psi[10])))
    abline(h=psi3[1], lty=2)
    title(main=expression(paste(psi[0], ' values')))
    
    boxplot(cbind(psi.acc.3.dtrreg[,1],
                  psi.acc.2.dtrreg[,1],
                  psi.acc.1.dtrreg[,1]),
            names=c(expression(psi[30]),
                    expression(psi[20]),
                    expression(psi[10])))
    abline(h=psi3[1], lty=2)
    title(main=expression(paste(psi[0], ' values using DTRreg')))
    
    png(paste('dwols-multi-stage-dtrreg-psi1_n', iter, '.png', sep=''), width=10, height=6, units='in', res=300)
    par(mfrow=c(1, 2))
    boxplot(cbind(psi.acc.3[,2],
                  psi.acc.2[,2],
                  psi.acc.1[,2]),
            names=c(expression(psi[31]),
                    expression(psi[21]),
                    expression(psi[11])))
    abline(h=psi3[2], lty=2)
    title(main=expression(paste(psi[1], ' values')))
    
    boxplot(cbind(psi.acc.3.dtrreg[,2],
                  psi.acc.2.dtrreg[,2],
                  psi.acc.1.dtrreg[,2]),
            names=c(expression(psi[31]),
                    expression(psi[21]),
                    expression(psi[11])))
    abline(h=psi3[1], lty=2)
    title(main=expression(paste(psi[1], ' values using DTRreg')))
    dev.off()
  }
}

# stop timer
proc.time() - ptm
