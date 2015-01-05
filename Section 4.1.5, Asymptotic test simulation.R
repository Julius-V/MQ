rm(list = ls())
source("func.R")
### Performing Monte Carlo simulation of 
### the asymptotic test for moving maximum model


library(foreach)
library(doSNOW)
cl <- makeCluster(2)
registerDoSNOW(cl)


# H_0
# Length of the auxiliary process (the one that gives parameter estimates)
N <- 15000
# Burn in
bi <- 1500
# Iterations
iter <- 1500

resultsH0 <- array(0, c(iter, 3))
n <- 300
M <- 15
thetas <- c(0.1, 0.5, 0.95)
for(thetaR in thetas) {
  resultsH0[, thetaR == thetas] <- 
    foreach(b = iter(1:iter), .combine = "c", .packages = "caTools") %dopar% {
      y <- mq.sim(n, M = M, pr = 1, theta = thetaR)
      asyTest(y, M = M, N = N, n.start = bi)
    }
}



# H_1
# Length of the auxiliary process (the one that gives parameter estimates)
N <- 10000
# Burn in
bi <- 1500
# Iterations
iter <- 1500

resultsH1 <- array(0, c(iter, 3))
n <- 250
M <- 15
phis <- c(0.1, 0.5, 0.95)
for(phi in phis) {
  resultsH1[, phi == phis] <- 
    foreach(b = iter(1:iter), .combine = "c", .packages = "caTools") %dopar% {
      y <- as.numeric(arima.sim(list(ar = phi), n))
      asyTest(y, M = M, N = N, n.start = 2000)
    }
}
stopCluster(cl)