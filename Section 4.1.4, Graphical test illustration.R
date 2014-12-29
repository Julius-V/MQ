rm(list = ls())
source("func.R")
### The code for the illustration of the graphical test

library(ggplot2)
library(np)
library(zoo)

## Simulation setting, four situations
# Parameter M
Ms <- c(20, 10, 20, 20)
# Sample size
Ns <- c(20000, 20000, 20000, 5000)
# Quantiles
ps <- c(1, 1, 0.7, 1)
# Bootstrap iterations
bit <- 250
# True parameter theta
thetaR <- 0.95


plotData <- list()
set.seed(123)
for(i in 1:4) {
  # Generating process with given parameters
  x <- mq.sim(Ns[i], 20, ps[i], const = 0, theta = thetaR, n.start = 5000, sd = 1)
  # Obtaining moving maximum
  qs <- turboQ(x, 1, Ms[i])[,1]
  # Regime lengths
  cnts <- tapply(qs, qs, length)[as.character(unique(qs))]
  # Levels
  vls <- as.numeric(names(cnts))
  # And maximal lengths
  maxLs <- sapply(vls, function(v)
    Ms[i] - (which.min(abs(qs - v)) - which.min(abs(tail(x, length(qs)) - v)))) + 1
  # Choosing only full length regimes
  idx <- maxLs == Ms[i]
  # Model
  mod <- mq(x, K = 0, M = Ms[i], p = 1, const = FALSE)
  theta <- mod$coef
  
  x2 <- cnts[idx]
  x1 <- vls[idx]
  x2 <- ordered(x2)
  # Estimating joint density
  dn <- npudens(~ x1 + x2)
  # A function to estimate parameter p at a given level
  npr <- function(x) {
    w <- predict(dn, newdata = data.frame(x1 = x, x2 = ordered(1:max(x2))))
    m <- sum(1:max(as.numeric(x2)) * w) / sum(w)
    momEst(m, Ms[i])
  }
  # Range of levels of interest
  from <- quantile(x1, 0.1)
  to <- quantile(x1, 0.9)
  
  # Bootstrap to get a confidence set
  bts1 <- replicate({
    # Sampling
    idx <- sample(length(x1), replace = TRUE)
    aux1 <- x1[idx]
    aux2 <- x2[idx]
    # Estimating a joint density with a given bandwidth (otherwise would take forever)
    dn0 <- npudens(~ aux1 + aux2, bandwidth.compute = FALSE, bws = dn$bw)
    # Again a function to estimate parameter p at a given level
    npr0 <- function(x) {
      w <- predict(dn0, newdata = data.frame(aux1 = x, aux2 = ordered(unique(as.numeric(x2)[idx]))))
      m <- sum(unique(as.numeric(x2)[idx]) * w) / sum(w)
      momEst(m, Ms[i])
    }
    # Computing it over the range of interest
    sapply(seq(from, to, length = 100), npr0)}, n = bit)
  
  # Obtaining a curve of p's
  zht <- sapply(seq(from, to, length = 100), npr)
  ## Part II - using model
  eps <- resid(mod)
  # Transforming the range of interest
  from <- from * (1-theta)
  to <- to * (1-theta)
  # Bootstrap to get a confidence set
  bts2 <- replicate({
    # Sampling
    sm <- sample(eps, replace = TRUE)
    # Obtaining p's    
    1 - ecdf(sm)(seq(from, to, length = 100))}, n = bit)  
  # Another curve of p's
  pn <- 1 - ecdf(eps)(seq(from, to, length = 100))
  
  ## Part III - results
  # Transforming the range of interest back to original units
  xs <- seq(1 / (1 - theta) * from, 1 / (1 - theta) * to, length = 100)
  plotData[[i]] <- data.frame(x = xs, y = c(zht, apply(bts1, 1, quantile, 0.025),
                                       apply(bts1, 1, quantile, 0.975), pn,
                                       apply(bts2, 1, quantile, 0.025),
                                       apply(bts2, 1, quantile, 0.975)),
                         confint = rep(c(0, 1, 1, 0, 1, 1), each = 100),
                         type = rep(c("Data", "Model"), each = 300), sep = rep(1:6, each = 100))
}
plotDf <- do.call(rbind, plotData)
plotDf <- cbind(plotDf, cs = rep(paste("Case", 1:4), each = 600))

ggplot(plotDf, aes(x = x, y = y, color = type, linetype = factor(confint), group = sep)) +
  geom_line() + theme_bw() + thm +
  scale_linetype(guide = "none") + scale_color_grey("Source", end = 0.6, start = 0.1) +
  xlab(expression(paste("Level ", x[t]))) + ylab("Structural break probability") + 
  facet_wrap(~cs, scales = "free") + theme(legend.position = 'bottom')