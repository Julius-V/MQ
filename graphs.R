rm(list = ls())
source("func.R")
library(ggplot2)

# Fixed M - Figure C.1
M <- 20
data <- data.frame(y = c(dtgeom(1:M, 0.01, M), dtgeom(1:M, 0.03, M), 
                         dtgeom(1:M, 0.05, M), dtgeom(1:M, 0.07, M)), 
                   x = 1:M, p = factor(rep(c(c(1, 3, 5, 7) / 100), each = M)))
ggplot(data, aes(y = y, x = x, group = p, color = p, linetype = p)) + 
  geom_line() + geom_point() + theme_bw() + ylab("P(X = k)") + xlab("k") +
  thm + scale_linetype("Parameter p") + scale_color_grey("Parameter p")


# Fixed p - Figure C.2
p <- 0.03
data <- data.frame(y = c(dtgeom(1:10, p, 10), dtgeom(1:20, p, 20), dtgeom(1:30, p, 30),
                         dtgeom(1:40, p, 40)), x = c(1:10, 1:20, 1:30, 1:40), 
                   M = factor(rep(c(10, 20, 30, 40), c(10, 20, 30, 40))))
ggplot(data, aes(y = y, x = x, group = M, color = M, linetype = M)) + 
  geom_line() + geom_point() + theme_bw() + ylab("P(X = k)") + xlab("k") +
  thm + scale_linetype("Parameter M") + scale_color_grey("Parameter M")


# Structural breaks example - Figure C.3
set.seed(123)
x <- mq.sim(5000, M = 50, theta = 0.95, pr = 1, n.start = 5000)
plotData <- data.frame(x = 1:5000, y = x)
ggplot(plotData, aes(x, y)) + geom_line() + theme_bw() + thm + ylab(NULL) + xlab("Time")


# Moving maximum processes, \xi* densities and their ACF's - Figures C.4, C.5
N <- 50000
Ms <- c(40, 40, 5, 5)
thetas <- c(0.95, 0.1, 0.95, 0.1)

set.seed(123)
plots <- sapply(1:4, function(n) {
  # Generating a process
  x <- mq.sim.aux(N, M = Ms[n], pr = 1, theta = thetas[n])
  # Obtaining moving quantiles
  qs <- turboQ(x[[1]], 1, Ms[n])[,1]
  # Process realization
  p1 <- ggplot(data.frame(y = x[[1]], x = 1:N), aes(x, y)) + geom_line() + thm + theme_bw() +
    xlab("Time") + ylab(expression("x"[t])) + 
    ggtitle(bquote(paste("Process realization, M=", .(Ms[n]), ", ", theta, "=", .(thetas[n]))))
  # Finding \xi*'s
  ns <- match(unique(qs), x[[1]])
  expr <- c(expression(paste('  ', xi[t])), expression(paste('  ', xi[t]^'*')))
  p2 <- ggplot(data.frame(x = c(x[[2]], x[[2]][ns]), id = factor(rep(1:2, c(length(x[[2]]), length(ns))))),
               aes(x, color = id, group = id, linetype = id)) + 
    geom_density() + thm + theme_bw() + xlab(NULL) + ylab(NULL) + 
    scale_linetype("", labels = expr) + scale_color_grey("", labels = expr, end = 0.5) +
    theme(legend.text = element_text(size = 15), legend.position = 'bottom') + ggtitle("Density functions")
  p3 <- qacf(x[[2]][ns], title = "Correlogram")
  list(p1, p2, p3)
})
library(gridExtra)
do.call(grid.arrange, c(plots[c(1,4,2,5,3,6)], list(ncol = 2)))
do.call(grid.arrange, c(plots[c(7,10,8,11,9,12)], list(ncol = 2)))


## Empirical ACF vs Approximated ACF - Figure C.6
set.seed(123)
# True parameters theta
thetas <- c(0.95, 0.95)
# Parameters M
Ms <- c(60, 20)
# 2 * Delta_1
sth <- c(14, 4)
# lag.max for acf()
nl <- 1000
plots <- list()
for(n in 1:2) {
  # Generating a process
  x <- mq.sim(100000, pr = 1, theta = thetas[n], M = Ms[n])
  # Obtaining moving quantiles
  qs <- turboQ(x, 1, Ms[n])[,1]
  # Extracting regime lengths
  cnts <- tapply(qs, qs, length)[as.character(unique(qs))]
  # And their levels
  vls <- unique(qs)
  # And their maximal lengths
  maxLs <- sapply(vls, function(v)
    Ms[n] - (which.min(abs(qs - v)) - which.min(abs(tail(x, length(qs)) - v)))) + 1
  # Empirical ACF
  bacf <- acf(x, plot = FALSE, lag.max = nl)
  # ...and the approximated one
  apprx <- thetas[n]^(sth[n] + 0:nl / (mean(cnts) + mean(Ms[n] - maxLs[maxLs <= Ms[n]])))
  bacfdf <- with(bacf, data.frame(lag, acf, another = apprx))
  p <- qplot(lag, acf, data = bacfdf, geom = "bar", stat = "identity",
             position = "identity", main = "") + theme_bw() +
    xlab(NULL) + geom_bar(fill = "grey50", stat = "identity") +
    thm + theme_bw() + ylab(NULL) + theme(legend.position = 'none') +
    geom_line(aes(x = lag, y = another), size = 1.2, color = "grey20")
  plots[[n]] <- p
}
library(gridExtra)
do.call(grid.arrange, c(plots, nrow = 1))