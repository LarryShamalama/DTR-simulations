####################
### G-estimation ###
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
#alpha <- c(-1, 3, 2)

num.iter <- c(100, 1000, 10000)
num.patients <- 1000

use.dtrreg <- FALSE # if analysis will be performed alongside DTRreg package

# Setup

create.data <- function(n){
  X <- rnorm(n=n)
  A <- rbinom(n, 1, 1/(1+exp(-X)))
  Y <- rnorm(n=n, mean=exp(X) - X^3 + A*(psi[1] + X*psi[2]), sd=1)
  
  return (cbind(Y, A, X))
}

psi.func <- function(.A, .Ahat, .H.beta, .H.psi, .Y){
  W <- diag(.A - .Ahat) - (.A - .Ahat)*.H.beta %*% solve(t(.H.beta) %*% .H.beta) %*% t(.H.beta)
  psi <- as.vector(solve( t(.H.psi) %*% W %*% (.A*.H.psi)) %*% t(.H.psi) %*% W %*% .Y)  
  
  return (psi)
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
    psi.model1   <- psi.func(A, Ahat.inc, x.beta.inc, x.psi, Y)
    
    # treatment model incorrect
    psi.model2   <- psi.func(A, Ahat.inc, x.beta.cor, x.psi, Y)
    
    # outcome model incorrect
    psi.model3   <- psi.func(A, Ahat.cor, x.beta.inc, x.psi, Y)
    
    # both models correct
    psi.model4   <- psi.func(A, Ahat.cor, x.beta.cor, x.psi, Y)
    
    return (rbind(psi.model1,
                  psi.model2,
                  psi.model3,
                  psi.model4))
  }
  
  ptm <- proc.time() # start timer
  obj <- foreach(x = 1:n, .combine=rbind) %dopar% f(x)
  
  indices <- seq(1, n)*4
  g.model1.psi <- obj[indices-3,]
  g.model2.psi <- obj[indices-2,]
  g.model3.psi <- obj[indices-1,]
  g.model4.psi <- obj[indices,]
  
  if (!use.dtrreg){
    png(paste('single-stage-psi_n', n, '.png', sep=''), width=10, height=6, units='in', res=300)
    par(mfrow=c(1, 2))
    boxplot(cbind(as.vector(g.model1.psi[,1]),
                  as.vector(g.model2.psi[,1]),
                  as.vector(g.model3.psi[,1]),
                  as.vector(g.model4.psi[,1])),
            ylim=c(-1, 0.8))
    abline(h=psi[1], lty=2, lwd=1)
    title(xlab='Model type',
          ylab=expression(paste(psi[0], '')),
          main=expression(paste(psi[0], ' values')))
    
    boxplot(cbind(g.model1.psi[,2],
                  g.model2.psi[,2],
                  g.model3.psi[,2],
                  g.model4.psi[,2]),
            ylim=c(-0.8, 4))
    abline(h=psi[2], lty=2, lwd=1)
    title(xlab='Model type',
          ylab=expression(paste(psi[1], '')),
          main=expression(paste(psi[1], ' values')))
    dev.off()
  }
  
  
  if (use.dtrreg){  
    
    dtr.model1.psi <- c()
    dtr.model2.psi <- c()
    dtr.model3.psi <- c()
    dtr.model4.psi <- c()

    f.dtrreg <- function(i){
      temp.data <- create.data(num.patients)
      Y <- temp.data[,1]
      A <- temp.data[,2]
      X <- temp.data[,3]
      
      X1 <- X^3
      X2 <- exp(X)
      
      blip.mod  <- list(~X) # has to be correct
      treat.mod.cor <- list(A~X)
      treat.mod.inc <- list(A~1)
      tf.mod.cor    <- list(~X1+X2)
      tf.mod.inc    <- list(~X)
      
      model1 <- DTRreg(Y, blip.mod, treat.mod.inc, tf.mod.inc, method='gest')
      model2 <- DTRreg(Y, blip.mod, treat.mod.inc, tf.mod.cor, method='gest')
      model3 <- DTRreg(Y, blip.mod, treat.mod.cor, tf.mod.inc, method='gest')
      model4 <- DTRreg(Y, blip.mod, treat.mod.cor, tf.mod.cor, method='gest')
      
      dtr.model1.psi <- rbind(dtr.model1.psi, model1$psi[[1]])
      dtr.model2.psi <- rbind(dtr.model2.psi, model2$psi[[1]])
      dtr.model3.psi <- rbind(dtr.model3.psi, model3$psi[[1]])
      dtr.model4.psi <- rbind(dtr.model4.psi, model4$psi[[1]])
    }
    
    ptm <- proc.time() # start timer
    obj <- foreach(x = 1:n, .combine=rbind) %dopar% f.dtrreg(x)
    
    indices <- seq(1, n)*4
    dtr.model1.psi <- obj[indices-3,]
    dtr.model2.psi <- obj[indices-2,]
    dtr.model3.psi <- obj[indices-1,]
    dtr.model4.psi <- obj[indices,]

    png(paste('single-stage-psi0_n', n, '.png', sep=''), width=10, height=6, units='in', res=300)
    par(mfrow=c(1, 2))
    boxplot(cbind(g.model1.psi[,1],
                  g.model2.psi[,1],
                  g.model3.psi[,1],
                  g.model4.psi[,1]),
            ylim=c(-1, 0.8))
    abline(h=psi[1], lty=2, lwd=1)
    title(xlab='Model type',
          ylab=expression(paste(psi[0], ' values')),
          main=paste('Single-stage\nn = ', n))
    boxplot(cbind(dtr.model1.psi[,1],
                  dtr.model2.psi[,1],
                  dtr.model3.psi[,1],
                  dtr.model4.psi[,1]),
            ylim=c(-1, 0.8))
    abline(h=psi[1], lty=2, lwd=1)
    title(xlab='Model type',
          ylab=expression(paste(psi[0], ' values')),
          main='Using DTRreg')
    dev.off()
    
    png(paste('single-stage-psi1_n', n, '.png', sep=''), width=10, height=6, units='in', res=300)
    par(mfrow=c(1, 2))
    boxplot(cbind(g.model1.psi[,2],
                  g.model2.psi[,2],
                  g.model3.psi[,2],
                  g.model4.psi[,2]),
            ylim=c(-0.8, 4))
    abline(h=psi[2], lty=2, lwd=1)
    title(xlab='Model type',
          ylab=expression(paste(psi[1], ' values')),
          main=paste('Single-stage\nn = ', n))
    boxplot(cbind(dtr.model1.psi[,2],
                  dtr.model2.psi[,2],
                  dtr.model3.psi[,2],
                  dtr.model4.psi[,2]),
            ylim=c(-0.8, 4))
    abline(h=psi[2], lty=2, lwd=1)
    title(xlab='Model type',
          ylab=expression(paste(psi[1], ' values')),
          main='Using DTRreg')
    dev.off()
  }
}

# stop timer
proc.time() - ptm
