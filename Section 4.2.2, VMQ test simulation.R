rm(list = ls())
source("func.R")
### Performing Monte Carlo simulation of 
### the VMQ specification test


library(foreach)
library(doSNOW)
cl <- makeCluster(2)
registerDoSNOW(cl)

# Iterations
it <- 2000
# Theoretical covariance matrix
sigma <- matrix(c(2, 0.6, 0.6, 1), byrow = TRUE, nrow = 2)
# Parameter M
M <- 20
coefs <- c(0.1, 0.2, 0.3)
# Quantiles
prs <- c(0.25, 1)
# Sample sizes
ns <- c(150, 300, 600, 1000)

results <- array(0, c(it, 4 * 2, 3))
resultsH1 <- array(0, c(it, 4, 2, 3))

# Under H_0 (very inefficient but does the job)
set.seed(123)
for(coeff in coefs) {
  phiM <- matrix(c(coeff, coeff / 2, coeff, coeff / 2,
                   coeff, coeff / 2, coeff, coeff / 2), nrow = 2, byrow = TRUE)
  results[, , coeff == coefs] <- 
    foreach(b = iter(1:it), .combine = 'rbind', .packages = "MASS") %dopar% {
      data <- lapply(ns, vmq.sim, M = M, pr = prs[1], phi = phiM, theta =  matrix(0, 2, 2), cov = sigma, n.start = 1000)
      vmods1 <- lapply(1:length(ns), function(n) vmq(data[[n]], M = M, K = 2, pr = prs[1], const = FALSE, cov = sigma))
      vmods2 <- lapply(1:length(ns), function(n) vmq(data[[n]], M = M, K = 2, pr = prs[2], const = FALSE, cov = sigma))
      c(sapply(vmods1, function(m) vmqTest(m)$stat), sapply(vmods2, function(m) vmqTest.vmq(m)$stat))
    }
}

# Under H_1
coeff <- 0.1
phiM <- matrix(c(coeff, coeff / 2, coeff, coeff / 2,
                 coeff, coeff / 2, coeff, coeff / 2), nrow = 2, byrow = TRUE)
for(pr in prs)
  for(coeff in coefs) {
    thetaM <- matrix(c(coeff, coeff, coeff, coeff), nrow = 2, byrow = TRUE)
    resultsH1[, , pr == prs, coeff == coefs] <- 
      foreach(b = iter(1:it), .combine = 'rbind', .packages = "MASS") %dopar% {
        data <- lapply(ns, vmq.sim, M = M, pr = pr, phi = phiM, theta =  thetaM, cov = sigma, n.start = 1000)
        vmods <- lapply(1:length(ns), function(n) vmq(data[[n]], M = M, K = 2, pr = pr, const = FALSE, cov = sigma))
        sapply(vmods, function(m) vmqTest.vmq(m)$stat)
      }
  }
stopCluster(cl)