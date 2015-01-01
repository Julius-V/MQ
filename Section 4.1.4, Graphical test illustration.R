rm(list = ls())
source("func.R")
### The code for the illustration of the graphical test
### Figure 4.1

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
bit <- 300
# True parameter theta
thetaR <- 0.95


plotData <- list()
set.seed(123)
for(i in 1:4) {
  print(i)
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
  # Model
  mod <- mq(x, K = 0, M = Ms[i], p = 1, const = FALSE)
  plotData[[i]] <- graphTest(cnts, vls, Ms[i], maxLs, boot.iter = bit, alpha = 0.05, model = mod)
}
plotDf <- do.call(rbind, plotData)
plotDf <- cbind(plotDf, cs = rep(paste("Case", 1:4), each = 600))

ggplot(plotDf, aes(x = x, y = y, color = type, linetype = factor(confint), group = sep)) +
  geom_line() + theme_bw() + thm +
  scale_linetype(guide = "none") + scale_color_grey("Source", end = 0.6, start = 0.1) +
  xlab(expression(paste("Level ", x[t]))) + ylab("Structural break probability") + 
  facet_wrap(~cs, scales = "free") + theme(legend.position = 'bottom')