rm(list = ls())
source("func.R")
### Simulation study to compare ML and MoM estimates

## Simulation setting
# Values of parameter p
ps <- c(0.01, 0.05, 0.1)
# Values of parameter M
Ms <- c(10, 30, 100)
# Samples sizes
ns <- c(10, 50, 150)


set.seed(123)
ar <- array(0, c(3, 4, 3, 3))
for(i in 1:3)
  for(j in 1:3)
    for(k in 1:3) {
      # "Defining" a particular truncated geometric r.v. X
      p <- ps[i]
      M <- Ms[j]
      n <- ns[k]
      dt <- replicate({
        ## Sampling from X
        # P(X = 1...M)=probs
        probs <- cumsum(c(p * (1 - p)^(0:(M - 2)), (1 - p)^(M - 1)))
        rns <- runif(n)
        # Sample
        smp <- findInterval(rns, probs) + 1
        # ML and MoM estimates
        c((n - sum(smp == M)) / (n - sum(smp == M) + n * (mean(smp) - 1)), momEst(mean(smp), M))
      }, n = 2000)
      ar[i, , j, k] <- c(rowMeans(dt), apply(dt, 1, sd))
    }


# Table C.1
tmp <- apply(ar, c(3, 1, 4), function(x) x[4] / x[3])
fin <- cbind(tmp[, , 1], tmp[, , 2], tmp[, , 3])
cat(paste(apply(round(fin, 2), 1, paste, collapse = " & "), collapse = "\n"))

# Table C.2
tmp <- apply(ar, c(3, 1, 4), function(x) {
  cn <- ps[which.min(abs(x[1] - ps))]
  abs(x[1] - cn) / abs(x[2] - cn)
})
fin <- cbind(tmp[, , 1], tmp[, , 2], tmp[, , 3])
cat(paste(apply(round(fin, 2), 1, paste, collapse = " & "), collapse = "\n"))