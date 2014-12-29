rm(list = ls())
source("func.R")
### An application using the moving maximum model on
### the absolute returns of the Bank of America

library(xtable)
library(np)
library(sandwich)
library(ggplot2)
library(gridExtra)
library(quantmod)

# Fixed parameter M
M <- 6
# Bootstrap iterations
bit  <- 400

# Data
x <- getSymbols("BAC", auto.assign = FALSE, from = "2011-10-17", to = "2014-07-30")[, 6]
# Transforming to absolute returns and setting intercept equal to zero (skipping computations)
tr <- function(x) abs(diff(log(x))) - 0.006868 / (1 - 0.256746)
y <- tr(x)[is.finite(tr(x))]
do.call(grid.arrange, list(qplot(x = index(y), y = as.numeric(y), geom = "line") + thm + theme_bw() +
                             xlab("Time") + ylab(NULL), qacf(as.numeric(y)), nrow = 1))

# Model
mod <- mq(as.numeric(y), pr = 1, M = M, K = 0, const = FALSE)
theta <- mod$coef
eps <- resid(mod)

do.call(grid.arrange, list(qplot(x = index(y[-1:-6]), y = eps, geom = "line") + thm + 
                             theme_bw() + xlab("Time") + ylab(NULL), qacf(eps), nrow = 1))

# Regimes analysis
y <- as.numeric(y)
# Moving maximum
qs <- turboQ(y, 1, M)[,1]
# Regime lengths
cnts <- tapply(qs, qs, length)[as.character(unique(qs))]
# Regime values
vls <- unique(qs)
# Maximal lengths
maxLs <- sapply(vls, function(v)
  M - (which.min(abs(qs - v)) - which.min(abs(tail(y, length(qs)) - v)))) + 1
# Estimates of p at each regime level
ps <- 1 - ecdf(eps)((1 - theta) * vls)
# Corresponding mean regime lengths
meanEst <- function(p, M) (2 * (1 - (1 - p)^M * (1 + p * M)) / p^2 - (1 - (1 - p)^M) / p)
means <- meanEst(ps, maxLs)
# Ignoring annomalies
idx <- maxLs >= 1 & maxLs <= M & cnts <= M


## Results
# Table 1
xtable(addmargins(table(cnts[idx], maxLs[idx])))

# Tests
mqTest(mod)
mqTest(mod, rob = TRUE)

# Another test
f <- var(cnts[idx]) * (1 + 2 * sum(acf(cnts[idx], plot = FALSE, lag.max = 2)$acf[-1]))
1/(sqrt(f)*sqrt(length(cnts[idx]))) * sum(cnts[idx] - means[idx])

### Graphical test - Figure 3
set.seed(123)
idx <- maxLs == M & cnts <= M
## Part I - using data only (no model)
x2 <- cnts[idx]
x1 <- vls[idx]
x2 <- ordered(x2)
# Estimating a joint density
dn <- npudens(~ x1 + x2)
# A function to estimate parameter p at a given level
npr <- function(x) {
  w <- predict(dn, newdata = data.frame(x1 = x, x2 = ordered(1:max(x2))))
  m <- sum(1:max(as.numeric(x2)) * w) / sum(w)
  momEst(m, M)
}
# Range of levels of interest
from <- quantile(x1, 0.1)
to <- quantile(x1, 0.9)

# Bootstrap to get a confidence set
iii <- 1
bts1 <- replicate({
  print(iii)
  iii <<- iii + 1
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
    momEst(m, M)
  }
  # Computing it over the range of interest
  sapply(seq(from, to, length = 100), npr0)}, n = bit)

# Obtaining a curve of p's
zht <- sapply(seq(from, to, length = 100), npr)
## Part II - using model
eps <- resid(mod)
# Transforming the range of interest
from <- from * (1 - theta)
to <- to * (1 - theta)
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
plotData <- data.frame(x = xs, y = c(zht, 
                                     apply(bts1, 1, quantile, 0.025),
                                     apply(bts1, 1, quantile, 0.975), 
                                     pn,
                                     apply(bts2, 1, quantile, 0.025),
                                     apply(bts2, 1, quantile, 0.975)),
                       confint = rep(c(0, 1, 1, 0, 1, 1), each = 100),
                       type = rep(c("Data", "Model"), each = 300), sep = rep(1:6, each = 100))

ggplot(plotData, aes(x = x, y = y, color = type, linetype = factor(confint), group = sep)) +
  geom_line() + theme_bw() + thm +
  scale_linetype(guide = "none") + scale_color_grey("Source", end = 0.6, start = 0.1) +
  xlab(expression(paste("Level ", x[t]))) + ylab("Structural break probability") + 
  theme(legend.position = 'bottom')

# Figure 4
library(plot3D)
mm <- mesh(1:M, seq(min(y), max(y), length = 30))
xx <- mm$x
yy <- mm$y
pn <- 1 - ecdf(eps)(yy)
z <- pn^(1 - 1 * (xx == M)) * (1 - pn)^(xx - 1)
library(lattice)
print(wireframe(z ~ xx + yy, data = data.frame(xx = xx, yy = yy, z = z), xlab = "Length", 
                ylab = expression(tilde("x")[t-1]), scale = list(arrows = FALSE), 
                zlab = "Pr.", screen = list(z = 140, x = -70),
                scpos = list(x = 9, y = 8), box = FALSE))