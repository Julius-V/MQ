rm(list = ls())
source("func.R")
### An application using VMQ model and the data on
### disposable income and consumption expenditure

library(ggplot2)
library(gridExtra)

# Data
data <- read.table("data/e1.dat", header = TRUE)
data <- apply(data, 2, function(x) diff(log(x)))


# Model
M <- 4
mod <- vmq(as.matrix(data)[, 2:3], K = 1, pr = 1, M = M)
# Table 5.2
coef(mod[[1]])
# Specification testing
vmqTest(mod)$p.val


## Graphs
dates <- seq(1960.25, 1982.75, by = 0.25)
# Income and consumption (Figures E.2, E.3)
do.call(grid.arrange, 
        list(qplot(x = dates, y = data[, 2], geom = "line") + 
               thm + theme_bw() + xlab("Time") + ylab(NULL), qacf(data[, 2]), nrow = 1))
do.call(grid.arrange, 
        list(qplot(x = dates, y = data[, 3], geom = "line") + 
               thm + theme_bw() + xlab("Time") + ylab(NULL), qacf(data[, 3]), nrow = 1))
# Residuals (Figures E.4, E.5)
do.call(grid.arrange, 
        list(qplot(x = dates[-1:-M], y = resid(mod[[1]])[,1], geom = "line") + 
               thm + theme_bw() + xlab("Time") + ylab(NULL), qacf(resid(mod[[1]])[,1]), nrow = 1))
do.call(grid.arrange, 
        list(qplot(x = dates[-1:-M], y = resid(mod[[1]])[,2], geom = "line") + 
               thm + theme_bw() + xlab("Time") + ylab(NULL), qacf(resid(mod[[1]])[,2]), nrow = 1))


# Covariance matrix
var(resid(mod[[1]])) * 10000


## Impulse response analysis
# VAR(4)
library(vars)
# n ahead
nah <- 15
varMOD <- VAR(data[, 2:3], p = M)
# VAR(M) irf
varShocks <- irf(varMOD, ortho = FALSE, boot = FALSE, n.ahead = nah)
nm <- c("Disposable income", "Consumption expenditure")

# Graphs
# Impulse to disposable income (Figure 5.4)
s1 <- data.frame(y = c(c(varShocks[[1]][[1]]), c(mq.irf(mod, shock = c(1, 0), n.ahead = nah)[-1:-3,])),
                 x = factor(0:nah), id = rep(c("VAR", "VMQ"), each = 2 * (nah + 1)), 
                 id2 = factor(rep(c(nm, nm), each = nah + 1)))
ggplot(s1, aes(x = x, y = y, group = id2, color = id2)) + geom_line() +
  thm + theme_bw() + facet_wrap(~ id) + xlab("") + theme(legend.position = 'bottom') + 
  scale_color_grey("Variable", start = 0.45, end = 0.75) + ylab(NULL)

# Impulse to consumption expenditure (Figure 5.5)
s2 <- data.frame(y = c(c(varShocks[[1]][[2]]), c(mq.irf(mod, shock = c(0, 1), n.ahead = nah)[-1:-3,])),
                 x = factor(0:nah), id = rep(c("VAR", "VMQ"), each = 2 * (nah + 1)), 
                 id2 = factor(rep(c(nm, nm), each = nah + 1)))
ggplot(s2, aes(x = x, y = y, group = id2, color = id2)) + geom_line() +
  thm + theme_bw() + facet_wrap(~ id) + xlab("") + theme(legend.position = 'bottom') + 
  scale_color_grey("Variable", start = 0.45, end = 0.75) + ylab(NULL)