rm(list = ls())
source("func.R")
### Performing Monte Carlo simulation of 
### the conditional-heteroskedasticity-robust test in small samples


library(foreach)
library(doSNOW)
cl <- makeCluster(2)
registerDoSNOW(cl)

# Burn-in
bi <- 1000
# Iterations
it <- 2000
# Parameter M
M <- 15
coeffs <- c(0.2, 0.45, 0.7)
# Sample sizes
ns <- c(150, 300, 600, 1000)


# Under H_0
resultsH0 <- array(0, c(it, 4 * 2, 3))
set.seed(123)
for(coef in coeffs) {
  resultsH0[, , coef == coeffs] <- 
    foreach(b = iter(1:it), .packages = c("caTools", "sandwich"), .combine = "rbind") %dopar% {
      x <- numeric(max(ns) + bi)
      norm <- rnorm(max(ns) + bi)
      x[1:M] <- norm[1:M]
      for(i in (M + 1):(bi + max(ns)))
        x[i] <- 0.1 * x[i - 1] + norm[i] * sqrt(0.05 + coef * x[i - 1]^2)
      x <- tail(x, -bi)
      out <- lapply(ns, function(n) 
        mq(x[1:n], M = M, pr = 1, K = 1, const = FALSE, HC = FALSE))
      c(sapply(out, function(x) mqTest.mq(x)$stat), sapply(out, function(x) mqTest(x, rob = TRUE)$stat))
    }
}

# Under H_1
resultsH1 <- array(0, c(it, 4 * 2, 3))
set.seed(123)
for(coef in coeffs) {
  resultsH1[, , coef == coeffs] <- 
    foreach(b = iter(1:it), .packages = c("caTools", "sandwich"), .combine = "rbind") %dopar% {
      x <- numeric(max(ns) + bi)
      norm <- rnorm(max(ns) + bi)
      x[1:M] <- norm[1:M]
      for(i in (M + 1):(bi + max(ns)))
        x[i] <- 0.1 * x[i - 1] + 0.1 * max(x[(i - 1):(i - M)]) + norm[i] * sqrt(0.05 + coef * x[i - 1]^2)
      x <- tail(x, -bi)
      out <- lapply(ns, function(n) 
        mq(x[1:n], M = M, pr = 1, K = 1, const = FALSE, HC = FALSE))
      c(sapply(out, function(x) mqTest.mq(x)$stat), sapply(out, function(x) mqTest(x, rob = TRUE)$stat))
    }
}
stopCluster(cl)