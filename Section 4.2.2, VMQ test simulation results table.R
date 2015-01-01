rm(list = ls())
load("data/vmqTestSim.RData")
### Constructing and printing tables from the MC simulation of the
### VMQ specification test
### Tables D.1 and D.2

alphas <- c(0.99, 0.95, 0.9)

# Results under H_0
tbl <- matrix(0, ncol = 3 * 3, nrow = 4 * 2)
for(n in 1:4)
  for(pr in 1:2)
    for(alpha in 1:3)
      for(coeff in 1:3)
        tbl[(n - 1) * 2 + pr, (coeff - 1) * 3 + alpha] <- 
          sprintf("%.2f", sum(results[, (pr - 1) * 4 + n, coeff] > qchisq(alphas[alpha], 2 * 2)) / dim(results)[1])
cat(paste(apply(tbl, 1, paste, collapse = " & "), collapse = "\n"))

# Results under H_1
for(n in 1:4)
  for(pr in 1:2)
    for(alpha in 1:3)
      for(coeff in 1:3)
        tbl[(n - 1) * 2 + pr, (coeff - 1) * 3 + alpha] <- 
          sprintf("%.2f", sum(resultsH1[, n, pr, coeff] > qchisq(alphas[alpha], 2 * 2)) / dim(resultsH1)[1])
cat(paste(apply(tbl, 1, paste, collapse = " & "), collapse = "\n"))