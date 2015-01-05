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
bit  <- 500

# Data
x <- getSymbols("BAC", auto.assign = FALSE, from = "2011-10-17", to = "2014-07-30")[, 6]
# Transforming to absolute returns and setting intercept equal to zero (skipping computations)
tr <- function(x) abs(diff(log(x))) - 0.006868 / (1 - 0.256746)
y <- tr(x)[is.finite(tr(x))]
# Figure 5.1
do.call(grid.arrange, list(qplot(x = index(y), y = as.numeric(y), geom = "line") + thm + theme_bw() +
                             xlab("Time") + ylab(NULL), qacf(as.numeric(y)), nrow = 1))

# Model
mod <- mq(as.numeric(y), pr = 1, M = M, K = 0, const = FALSE)
theta <- mod$coef
eps <- resid(mod)

# Figure E.1
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
# Ignoring annomalies
idx <- maxLs >= 1 & maxLs <= M & cnts <= M


## Results
# Table 5.1
xtable(addmargins(table(cnts[idx], maxLs[idx])))

# Tests
mqTest(mod)
mqTest(mod, rob = TRUE)

# Asymptotic test
library(fGarch)
(pars <- sstdFit(eps)$estim)
distr <- function(n) rsstd(n, pars[1], pars[2], pars[3], pars[4])
set.seed(123)
asyTest(y, M, distr = distr, N = 5e5)
# [1] 0.4503901
plotData <- data.frame(x = c(eps, distr(2e6)),
                       id = rep(c("Residuals", "Skew normal"), c(length(eps), 2e6)))
ggplot(plotData, aes(x = x, fill = id)) + geom_density(alpha = 0.15) + 
  theme_bw() + thm + xlab(NULL) + ylab(NULL) + xlim(c(min(eps), max(eps))) +
  theme(legend.position = "bottom") + scale_fill_grey("")


### Graphical test - Figure 5.2
set.seed(123)
plotData <- graphTest(cnts, vls, M, maxLs, boot.iter = bit, alpha = 0.05, model = mod)

ggplot(plotData, aes(x = x, y = y, color = type, linetype = factor(confint), group = sep)) +
  geom_line() + theme_bw() + thm +
  scale_linetype(guide = "none") + scale_color_grey("Source", end = 0.6, start = 0.1) +
  xlab(expression(paste("Level ", x[t]))) + ylab("Structural break probability") + 
  theme(legend.position = 'bottom')

# Figure 5.3
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