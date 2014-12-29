rm(list = ls())
load("data/robustTestSim.RData")
### Constructing and printing tables from the MC simulation of the
### conditional-heteroskedasticity-robust test for the MQ model

alphas <- c(0.99, 0.95, 0.9)

getEntry <- function(m, test, coeff, alpha)
  sprintf("%.2f", sum(m[, (test - 1) * 4 + n, coeff] > qchisq(alphas[alpha], 1)) / dim(m)[1])

# Results under H_0
tbl <- matrix(0, ncol = 3 * 3, nrow = 4 * 2)
for(n in 1:4)
  for(test in 1:2)
    for(alpha in 1:3)
      for(coeff in 1:3)
        tbl[(n - 1) * 2 + test, (coeff - 1) * 3 + alpha] <- getEntry(resultsH0, test, coeff, alpha)
cat(paste(apply(tbl, 1, paste, collapse = " & "), collapse = "\n"))


# Results under H_1
tbl <- matrix(0, ncol = 3 * 3, nrow = 4 * 2)
for(n in 1:4)
  for(test in 1:2)
    for(alpha in 1:3)
      for(coeff in 1:3)
        tbl[(n - 1) * 2 + test, (coeff - 1) * 3 + alpha] <- getEntry(resultsH1, test, coeff, alpha)
cat(paste(apply(tbl, 1, paste, collapse = " & "), collapse = "\n"))