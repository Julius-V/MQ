library(MASS)
library(caTools)

#' Random sampling of a truncated random variable
#' 
#' Random sampling of a truncated random variable with parameters `p` and `M`
#' @param n number of observations
#' @param p probabiliy of success
#' @param M positive integer, maximal possible length, truncation parameter
#' @return a vector of random draws of length `n`
#' @author Julius Vainora
rtgeom <- function(n, p, M) {
  pr <- c(p * (1 - p)^(0:(M - 2)), (1 - p)^(M - 1))
  rdraws <- runif(n)
  findInterval(rdraws, cumsum(c(0, pr)))
}

#' Probability mass function of a trunctated random variable
#' 
#' Probability mass function of a trunctated random variable with parameters `p` and `M`
#' 
#' @param x vector of quantities
#' @param p probability of success
#' @param M positive integer, maximal possible length, truncation parameter
#' @return a vector of probability mass values
#' @author Julius Vainora
dtgeom <- function(x, p, M) {
  pr <- c(p * (1 - p)^(0:(M - 2)), (1 - p)^(M - 1))
  pr[x]
}

#' Sligthly faster version of findInterval
mi <- function(x, vec)
  .Internal(findInterval(vec, x, FALSE, FALSE))

#' A function to compute a quantile when everything else is given
getQ <- function(s, h, hf = floor(h)) s[hf] + (h - hf) * (s[pmin(h + 1, length(s))] - s[hf])

#' Moving quantiles
#' 
#' Computing moving quantiles of probability `p` and window `M`
#' @param x numeric vector whose moving quantiles are wanted
#' @param p a vector of probabilities between 0 and 1
#' @param M window size for moving quantiles to compute
#' @return a matrix of `length(p)` columns and `length(x)-M` rows
#' @author Julius Vainora
turboQ <- function(x, p, M) {
  N <- length(x)
  h <- (M - 1) * p + 1
  s <- sort(x[1:M])
  qs <- matrix(NA, ncol = length(p), nrow = length(x) - M)
  for(i in (M + 1):length(x)) {
    qs[i - M, ] <- getQ(s, h)
    s <- s[-which(s == x[i - M])[1]]
    np <- mi(x[i], s)
    s <- c(s[seq_len(np)], x[i], s[setdiff(1:length(s), seq_len(np))])
  }
  qs
}

##' Fitting Moving Quantile Models
##'
##' `mq` estimates Moving Quantile regression using OLS and is a wrapper for \code{lm}.
##'
##' @param x numeric vector.
##' @param K nonnegative integer, the order of autoregressive term of the model.
##' @param M positive integer, window size used to compute quantiles.
##' @param pr numeric vector of probabilities with values in \emph{[0,1]}.
##' @param const logical, indicates whether constant term should be included.
##' @param type integer between 1 and 9 selecting one of the nine quantile algorithms detailed in \code{\link{quantile}}.
##' @param ... further arguments passed to \code{\link{lm}} function.
##' @return \code{mq} and/or \code{\link{lm}} object.
##' @author Julius Vainora
##' 
mq <- function(x, K = 0, M, pr = 0:4 / 4, const = TRUE, type = 7, ql = 0, HC = FALSE, vcovType = "HC3", ...) {
  # nera patikrinimu del ql
  if(!is.vector(x) | !is.numeric(x))
    stop("'x' is not a numeric vector")
  N <- length(x)
  L <- length(pr)
  if(N < L)
    stop("There have to be more observations than probabilities specified")
  if(M < L)
    stop("Window length has to be at least as long as 'pr'")
  if(M > N)
    stop("Too large window length")
  if(missing(K) | !is.numeric(K) | K < 0) {
    warning("'K' has to be nonnegative integer")
    K <- 0
  }
  if(missing(M) | !is.numeric(M) | M < 1)
    stop("'M' has to be positive integer") 
  if(!is.logical(const))
    stop("'const' should be logical")
  l <- max(K, M)
  X <- embed(x, l)[-(N - l + 1), , drop = FALSE]
  if(K > 0) {
    PHI <- X[, 1:K, drop = FALSE]
    colnames(PHI) <- paste0("ar", 1:K)
  } else 
    PHI <- NULL
  xm <- x[(l + 1):N]
  dt <- data.frame(cbind(x = xm, const = if(const) 1, PHI))
  if(is.null(pr))
    return(lm(x ~ -1 + ., dt))
  Q <- turboQ(x, pr, M)
  Q <- runquantile(x, M, pr, endrule = "trim")
  if(is.na(ncol(Q)))
    Q <- matrix(Q[-length(Q)])
  else
    Q <- Q[-nrow(Q), ]
  dim(Q) <- range(dim(Q))[2:1]
  if(any(ql != 0))
    Q <- sapply(1:length(ql), function(n) if(ql[n] != 0) c(rep(NA, ql[n]), Q[-(1:ql[n]), n]) else Q[, n])
  colnames(Q) <- paste0("q_", pr, "_", ql)
  out <- lm(x ~ -1 + ., data.frame(dt, Q), ...)
  oHC <- if(HC) vcovHC(out, type = vcovType) else NULL
  out <- c(out, x = list(xm), K = K, M = M, pr = list(pr), 
           Q = list(Q), PHI = list(cbind(if(const) rep(1, length(xm)), PHI)), HC = list(oHC))
  class(out) <- c("mq", "lm")
  out
}

##' A test for the MQ terms
##'
##' A test for the MQ terms in a given model.
##'
##' @param x an object of class "mq".
##' @param ... further arguments to be passed to or from methods.
##' @return \code{\link{lm}} object.
##' @author Julius Vainora
##' 
mqTest <- function(x, ...) {
  UseMethod("mqTest")
}

mqTest.default <- function(x, K = 0, M, pr = 0:4 / 4, const = TRUE, type = 7, ...)
  mqTest(mq(x, K, M, pr, const, type, ...))

mqTest.mq <- function(m, rob = FALSE) {
  PHI <- m$PHI
  if(!is.null(PHI))
    M <- diag(nrow(m$Q)) - PHI %*% solve(crossprod(PHI)) %*% t(PHI)
  else {
    M <- diag(nrow(m$Q))
    PHI <- rep(0, length(m$x))
  }
  H <- M %*% m$Q
  mod <- lm(m$x ~ -1 + PHI)
  r <- resid(mod)
  mod <- lm(r ~ -1 + H)
  theta <- coef(mod)
  if(rob)
    stat <- t(theta) %*% solve(vcovHC(mod)) %*% theta
  else
    stat <- c(crossprod(H %*% theta) / var(r))
  df <- length(theta)
  pval <- 1 - pchisq(stat, length(theta))
  out <- list(statistic = stat, parameter = df, 
              p.value = pval, theta = theta)
  class(out) <- "mqtest"
  out
}


##' Simulate from a Moving Quantile model
##'
##' Simulate from a Moving Quantile model
##'
##' @param n positive integer, length of output series.
##' @param M positive integer, window size used to compute quantiles.
##' @param pr numeric vector of probabilities with values in \emph{[0,1]}.
##' @param const numeric, constant term.
##' @param phi numeric vector, coefficients of autoregressive terms
##' @param theta numeric vector, coefficients of moving quantiles
##' @param type integer between 1 and 9 selecting one of the nine quantile algorithms detailed in \code{\link{quantile}}.
##' @param rand.gen optional, a function to generate the innovations.
##' @param innov an optional numeric vector of innovations. If not provided, \code{rand.gen} is used.
##' @param n.start positive integer, length of 'burn-in' period.
##' @param start.innov an optional numeric vector of innovations to be used for the burn-in period
##' @param ... additional arguments for \code{rand.gen}. Most usefully, the standard deviation of the innovations generated by rnorm can be specified by sd.
##' @return \code{\link{lm}} object.
##' @author Julius Vainora
##' 
mq.sim <- function(n, M, pr, const = 0, phi = 0, theta, type = 7, 
                   rand.gen = rnorm, innov = rand.gen(n, ...), n.start = 5000, 
                   start.innov = rand.gen(n.start, ...), ...) {
  if(missing(n) | !is.numeric(n) | length(n) > 1 || n < 1)
    stop("'n' has to be positive integer")
  if(missing(M) | !is.numeric(M) | length(M) > 1 || M < 1)
    stop("'M' has to be positive integer") 
  N <- n.start + n
  if(M > N)
    stop("Too large window size")
  if(!is.numeric(const) | length(const) > 1)
    stop("'const' has to be numeric")
  if(!is.numeric(phi))
    stop("'phi' has to be numeric vector")
  if(missing(theta))
    stop("Unspecified coefficients 'theta'")
  if(!is.numeric(theta))
    stop("'theta' has to be numeric vector")
  if(!is.numeric(n.start) | length(n.start) > 1 || n.start < 0)
    stop("'n.start' has to be nonnegative integer")
  out <- numeric(N)
  out[1:M] <- start.innov[1:M]
  tinnov <- c(start.innov, innov)
  for(i in (M + 1):N)
    out[i] <- const + phi %*% out[(i - length(phi)):(i - 1)] + 
    theta %*% quantile(out[(i - M):(i - 1)], pr, type = type) + tinnov[i]
  out[-1:-n.start]
}

graphTest <- function(cnts, vls, M, maxLs, lower = 0.1, upper = 0.9, n.points = 100, boot.iter = 300, alpha, model) {
  # Choosing only full length regimes and no annomalies
  idx <- maxLs == M & cnts <= M & cnts >= 1
  
  x2 <- cnts[idx]
  x1 <- vls[idx]
  x2 <- ordered(x2)
  # Estimating joint density
  dn <- npudens(~ x1 + x2)
  # A function to estimate parameter p at a given level
  npr <- function(x) {
    w <- predict(dn, newdata = data.frame(x1 = x, x2 = ordered(1:max(x2))))
    m <- sum(1:max(as.numeric(x2)) * w) / sum(w)
    momEst(m, M)
  }
  # Range of levels of interest
  from <- quantile(x1, lower)
  to <- quantile(x1, upper)
  
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
      momEst(m, M)
    }
    # Computing it over the range of interest
    sapply(seq(from, to, length = n.points), npr0)}, n = boot.iter)
  
  # Obtaining a curve of p's
  zht <- sapply(seq(from, to, length = n.points), npr)
  ## Part II - using model
  eps <- resid(model)
  # Transforming the range of interest
  theta <- model$coef
  from <- from * (1 - theta)
  to <- to * (1 - theta)
  # Bootstrap to get a confidence set
  bts2 <- replicate({
    # Sampling
    sm <- sample(eps, replace = TRUE)
    # Obtaining p's    
    1 - ecdf(sm)(seq(from, to, length = n.points))}, n = boot.iter)  
  # Another curve of p's
  pn <- 1 - ecdf(eps)(seq(from, to, length = n.points))
  
  ## Part III - results
  # Transforming the range of interest back to original units
  xs <- seq(1 / (1 - theta) * from, 1 / (1 - theta) * to, length = n.points)
  data.frame(x = xs, y = c(zht, 
                           apply(bts1, 1, quantile, alpha / 2),
                           apply(bts1, 1, quantile, 1 - alpha / 2), 
                           pn,
                           apply(bts2, 1, quantile, alpha / 2),
                           apply(bts2, 1, quantile, 1 - alpha / 2)),
             confint = rep(c(0, 1, 1, 0, 1, 1), each = n.points),
             type = rep(c("Data", "Model"), each = n.points * 3), 
             sep = rep(1:6, each = n.points))
}


HQC <- function(m) {
  k <- length(m$coeff)
  res <- resid(m)
  RSS <- sum(res^2)
  n <- length(res)
  n * log(RSS / n) + 2 * k * log(log(n))
}


##' MQ order selection
##'
##' The function returns information criteria for MQ models with specified AR order, window width and probabilies.
##'
##' @param x numeric vector.
##' @param K nonnegative integer, the order of autoregressive term of the model.
##' @param M positive integer, window size used to compute quantiles.
##' @param pr numeric vector of probabilities with values in \emph{[0,1]} to calculate all posible combinations from.
##' @param pr.fix list of fixed probability combinations.
##' @param nbest positive integer, the number of best combinations for each information criterion to return
##' @param IC characte vector of information criteria to be used
##' @param ... further arguments passed to \code{\link{mq}} function.
##' @return \code{\link{lm}} object.
##' @author Julius Vainora
##' 
mq.select <- function(x, K = 0, M, pr, pr.fix = NULL, prl = NULL,
                      nbest = 3, IC = c("AIC", "BIC", "HQC"), ...) {
  if(missing(pr) & missing(pr.fix))
    stop("Either 'pr' or 'pr.fix' must be supplied")
  if(!is.numeric(nbest) || nbest < 1)
    stop("'nbest' has to be a positive integer")
  IC <- intersect(IC, c("AIC", "BIC", "HQC"))
  matToList <- function(m)
    lapply(seq_len(ncol(m)), function(i) m[, i])
  if(!missing(pr)) {
    cmb <- lapply(lapply(1:length(pr), function(w) combn(pr, w)), matToList)
    cmb <- c(unlist(cmb, recursive = FALSE), list(NULL))
  } else cmb <- pr.fix
  names(prl) <- pr
  fl <- do.call(expand.grid, c(list(probs = 1:length(cmb), M = 1:length(M), K = 1:length(K)), prl))
  m <- t(apply(fl, 1, function(w) {
    if(M[w[2]] >= length(cmb[[w[1]]])) {
      mod <- mq(x, K = K[w[3]], M = M[w[2]], pr = cmb[[w[1]]], ql = w[4:length(w)][pr %in% w[1]])
      sapply(IC, function(f) get(f)(mod))
    } else rep(NA, length(IC))
  }))
  if(nrow(m) == 1) m <- t(m)
  fl <- data.frame(probs = sapply(cmb, paste, collapse = ", "), 
                   M = M[fl[, 2]], K = K[fl[, 3]], prl = fl[, 4:ncol(fl)])
  for(cr in 1:length(IC)){
    cat(IC[cr], "\n")
    idx <- order(m[, cr])[1:nbest]
    print(data.frame(fl[idx, ], value = m[idx, cr]), row.names = FALSE)
  }
  invisible(list(fl, m))
}

summary.mq <- function(x, ...) {
  sm <- 
    cat("Sum of coefficients:", sum(abs(coef(x))), "\n")
  summary.lm(x)
}

plot.mq <- function(m, ids = 1:ncol(m$Q), ...) {
  Q <- m$Q
  n <- floor(length(x) / 5)
  plot(m$x, type = 'l', ylab = "x")
  for(i in ids)
    lines(Q[, i], col = i + 1)
  par(ask = TRUE)
  acf(m$x, lag.max = n)
  acf(resid(m), lag.max = n)
  par(ask = FALSE)
}



##' Simulate from a Vector Moving Quantile model
##'
##' Simulate from a Vector Moving Quantile model
##'
##' @param n positive integer, length of output series.
##' @param M positive integer, window size used to compute quantiles.
##' @param pr numeric vector of probabilities with values in \emph{[0,1]}.
##' @param const numeric vector, constant terms.
##' @param phi numeric matrix, coefficients of autoregressive terms
##' @param theta numeric matrix, coefficients of moving quantiles
##' @param type integer between 1 and 9 selecting one of the nine quantile algorithms detailed in \code{\link{quantile}}.
##' @param rand.gen optional, a function to generate the innovations.
##' @param innov an optional numeric matrix of innovations. If not provided, \code{rand.gen} is used.
##' @param n.start positive integer, length of 'burn-in' period.
##' @param start.innov an optional numeric matrix of innovations to be used for the burn-in period
##' @param ... additional arguments for \code{rand.gen}. Most usefully, the standard deviation of the innovations generated by rnorm can be specified by sd.
##' @return \code{\link{lm}} object.
##' @author Julius Vainora
##' 
vmq.sim <- function(n, M, pr, const = 0, phi = matrix(0, nrow = nrow(theta)), theta, type = 7, 
                    cov = diag(nrow(theta)),
                    rand.gen = mvrnorm, innov = rand.gen(n, mu = rep(0, nrow(theta)), Sigma = cov, ...),
                    n.start = 5000, start.innov = rand.gen(n.start, mu = rep(0, nrow(theta)), Sigma = cov, ...), ...) {
  if(missing(n) | !is.numeric(n) | length(n) > 1 || n < 1)
    stop("'n' has to be positive integer")
  if(missing(M) | !is.numeric(M) | length(M) > 1 || M < 1)
    stop("'M' has to be positive integer") 
  N <- n.start + n
  if(M > N)
    stop("Too large window size")
  if(!is.numeric(const))
    stop("'const' has to be numeric vector")
  if(!is.matrix(phi) | !is.numeric(phi))
    stop("'phi' has to be numeric matrix")
  if(missing(theta))
    stop("Unspecified coefficients 'theta'")
  if(!is.matrix(theta) || !is.numeric(theta))
    stop("'theta' has to be numeric matrix")
  if(!is.numeric(n.start) | length(n.start) > 1 || n.start < 0)
    stop("'n.start' has to be nonnegative integer")
  out <- matrix(NA, nrow = N, ncol = nrow(theta))
  out[1:M, ] <- start.innov[1:M, ]
  tinnov <- rbind(start.innov, innov)
  for(i in (M + 1):N)
    out[i, ] <- const + phi %*% matrix(out[(i - 1):(i - ncol(phi) / nrow(phi)), ]) + 
    theta %*% matrix(apply(out[(i - M):(i - 1), ], 2, quantile, pr, type = type)) + tinnov[i, ]
  out[-1:-n.start, ]
}





##' Fitting Vector Moving Quantile Models
##'
##' `mq` estimates Vector Moving Quantile regression using OLS and is a wrapper for \code{lm}.
##'
##' @param x numeric matrix
##' @param K nonnegative integer, the order of autoregressive term of the model.
##' @param M positive integer, window size used to compute quantiles.
##' @param pr numeric vector of probabilities with values in \emph{[0,1]}.
##' @param const logical, indicates whether constant term should be included.
##' @param type integer between 1 and 9 selecting one of the nine quantile algorithms detailed in \code{\link{quantile}}.
##' @param ... further arguments passed to \code{\link{lm}} function.
##' @return \code{mq} and/or \code{\link{lm}} object.
##' @author Julius Vainora
##' 
vmq <- function(x, K = 0, M, pr = 0:4 / 4, const = TRUE, type = 7, cov = NULL, ...) {
  if(!is.matrix(x) | !is.numeric(x))
    stop("'x' is not a numeric matrix")
  N <- nrow(x)
  R <- ncol(x)
  L <- length(pr)
  if(N < L)
    stop("There have to be more observations than probabilities specified")
  if(M < L)
    stop("Window length has to be at least as long as 'pr'")
  if(M > N)
    stop("Too large window length")
  if(missing(K) | !is.numeric(K) | K < 0) {
    warning("'K' has to be nonnegative integer")
    K <- 0
  }
  if(missing(M) | !is.numeric(M) | M < 1)
    stop("'M' has to be positive integer") 
  if(!is.logical(const))
    stop("'const' should be logical")
  #   browser()
  l <- max(K, M)
  X <- embed(x, l)[-(N - l + 1), , drop = FALSE]
  if(K > 0) {
    PHI <- X[, 1:(R * K), drop = FALSE]
    colnames(PHI) <- paste0(" AR(", 1:R, ", t-", rep(1:K, each = R), ")")
    Z <- PHI
  } else 
    Z <- NULL
  xm <- x[(l + 1):N, ]
  if(const)
    Z <- cbind(rep(1, nrow(xm)), Z)
  Z <- t(Z)
  if(is.null(pr)) {
    mod <- lm(xm ~ -1 + t(Z))
    out <- c(model = list(mod), x = list(xm), K = K, M = M, pr = list(pr), cov = list(cov),
             const = if(const) list(matrix(coef(mod)[1, ])) else list(matrix(0, nrow = R)), type = type, Z = list(Z),
             theta = list(t(tail(coef(mod), R * L))), phi = list(t(coef(mod)[1 * const + 1:(R * K), ])))
    class(out) <- c("vmq", "lm")
    return(out)
  }
  Q <- t(lapply(1:R, function(r) {
    apply(X[, c(rep(FALSE, r - 1), TRUE, rep(FALSE, R - r)), drop = FALSE], 1, quantile, pr, type = type)
  }))
  Q <- do.call(rbind, Q)
  Q <- t(Q)
  colnames(Q) <- paste0(" q(p=", pr, ", ", rep(1:R, each = length(pr)), ")")
  data <- cbind(t(Z), Q)
  mod <- lm(xm ~ -1 + data)
  out <- c(model = list(mod), x = list(xm), K = K, M = M, pr = list(pr), cov = list(cov),
           const = if(const) list(matrix(coef(mod)[1, ])) else list(NULL), type = type, Q = list(Q), Z = list(Z),
           theta = list(t(tail(coef(mod), R * L))), phi = list(t(coef(mod)[1 * const + 1:(R * K), ])))
  class(out) <- c("vmq", "lm")
  out
}

library(ggplot2)
thm <-  theme(panel.grid.major = element_line(color = "grey95"),
              panel.grid.minor = element_line(color = "grey98"),
              legend.text = element_text(size = 11),
              legend.title = element_text(size = 11),
              axis.title.y = element_text(size = 12))

mq.irf <- function(x, ortho = TRUE, n.ahead = 10, cumulative = FALSE, shock) {
  pr <- x$pr
  M <- x$M
  theta <- x$theta
  phi <- x$phi
  const <- x$const
  R <- nrow(phi)
  K <- ncol(phi) / R
  L <- ncol(theta) / R
  l <- max(K, M)
  nm <- order(paste(1:R, rep(1:K, each = R)))
  
  PHI <- lapply(1:R, function(r) {
    phii <- c(sapply(split(phi[r, nm], rep(1:R, each = K)), function(w) c(w, rep(0, l - length(w)))))
    deltai <- matrix(0, nrow = l - 1, ncol = R * l)
    deltai[, (l * (r - 1) + 1):(l * r - 1)] <- diag(l - 1)
    rbind(phii, deltai)
  })
  PHI <- do.call(rbind, PHI)
  THETA <- lapply(1:R, function(r) rbind(theta[r, ], matrix(0, nrow = l - 1, ncol = R * L)))
  THETA <- do.call(rbind, THETA)
  if(length(THETA) == 0) {
    THETA <- 0
    impt <- 0
  } else 
    impt <- t(matrix(rapply(apply(theta, 1, split, rep(1:R, each = L)), sum), R, R))
  impp <- t(matrix(rapply(apply(phi[, nm], 1, split, rep(1:R, each = K)), sum), R, R))
  equa <- solve(diag(R) - impt - impp, const)
  out <- matrix(rep(equa, each = (l + n.ahead) * l), nrow = l + n.ahead, ncol = R * l)
  out[l, 1 + 0:(R - 1) * l] <- shock + out[l, 1 + 0:(R - 1) * l]
  cvec <- numeric(R * l)
  cvec[1 + 0:(R - 1) * l] <- const
  if(!is.null(pr)) {
    for(i in (l + 1):nrow(out))
      out[i, ] <- cvec + PHI %*% c(out[(i - 1):(i - l), 1 + 0:(R - 1) * l]) + 
      THETA %*% c(sapply(1:R, function(r) quantile(out[(i - 1):(i - M), 1 + (r - 1) * l], pr)))
  } else
    for(i in (l + 1):nrow(out))
      out[i, ] <- cvec + PHI %*% c(out[(i - 1):(i - l), 1 + 0:(R - 1) * l])
  if(cumulative)
    apply(out[, 1 + 0:(R - 1) * l], 2, cumsum)
  else
    out[, 1 + 0:(R - 1) * l]
}

vmqTest <- function(x, ...) {
  UseMethod("vmqTest")
}

vmqTest.vmq <- function(m) {
  Z <- t(m$Z)
  M <- diag(nrow(Z)) - Z %*% solve(crossprod(Z)) %*% t(Z)
  H <- M %*% m$Q
  r <- resid(lm(m$x ~ -1 + Z))
  R <- ncol(r)
  if(is.null(m$cov))
    m$cov <- var(r)
  theta <- c(apply(r, 2, function(x) solve(crossprod(H)) %*% t(H) %*% x))
  stat <- t(theta) %*% kronecker(solve(m$cov), t(H) %*% H) %*% theta
  df <- length(theta)
  pval <- 1 - pchisq(stat, length(theta))
  out <- list(statistic = stat, parameter = df, 
              p.value = pval, theta = theta)
  class(out) <- "vmqtest"
  out
}


mq.sim.aux <- function(n, M, pr, const = 0, phi = 0, theta, type = 7, 
                       rand.gen = rnorm, innov = rand.gen(n, ...), n.start = 5000, 
                       start.innov = rand.gen(n.start, ...), ...) {
  if(missing(n) | !is.numeric(n) | length(n) > 1 || n < 1)
    stop("'n' has to be positive integer")
  if(missing(M) | !is.numeric(M) | length(M) > 1 || M < 1)
    stop("'M' has to be positive integer") 
  N <- n.start + n
  if(M > N)
    stop("Too large window size")
  if(!is.numeric(const) | length(const) > 1)
    stop("'const' has to be numeric")
  if(!is.numeric(phi))
    stop("'phi' has to be numeric vector")
  if(missing(theta))
    stop("Unspecified coefficients 'theta'")
  if(!is.numeric(theta))
    stop("'theta' has to be numeric vector")
  if(!is.numeric(n.start) | length(n.start) > 1 || n.start < 0)
    stop("'n.start' has to be nonnegative integer")
  out <- numeric(N)
  out[1:M] <- start.innov[1:M]
  tinnov <- c(start.innov, innov)
  for(i in (M + 1):N)
    out[i] <- const + phi %*% out[(i - length(phi)):(i - 1)] + 
    theta %*% quantile(out[(i - M):(i - 1)], pr, type = type) + tinnov[i]
  list(out[-1:-n.start], innov)
}

##' Method of moments estimator
##' 
##' `momEst` numerically solves for a Method of Moments (MoM) estimate of a truncated geometric random variable parameter p
##' 
##' @param mu numeric, empirical mean of a truncated geometric random variable or regime lengths
##' @param M positive integer, maximal possible regime length
##' @return a length-one vector containing the MoM estimate
##' @author Julius Vainora
momEst <- function(mu, M)
  optim(0.5, function(p) ((1 - (1 - p)^M) / p - mu)^2, method = "BFGS")$par


qacf <- function(x, y = NULL, conf.level = 0.95, title = "") {
  ciline <- qnorm((1 - conf.level) / 2) / sqrt(length(x))
  bacf <- acf(x, plot = FALSE)
  bacfdf <- with(bacf, data.frame(lag, acf))
  significant <- (abs(bacfdf[, 2]) > abs(ciline))^2
  bacfdf <- cbind(bacfdf, significant)
  p <- qplot(lag, acf, data = bacfdf, geom = "bar", stat = "identity",
             position = "identity", main = title,
             fill = factor(-significant)) + theme_bw() +
    geom_hline(yintercept = -ciline, color = "black", size = 0.2, linetype = "dashed") +
    geom_hline(yintercept = ciline, color = "black", size = 0.2, linetype = "dashed") + 
    xlab(NULL) + geom_hline(yintercept = 0, color = "grey", size = 0.3) + 
    thm + theme_bw() + ylab(NULL) + scale_fill_grey("") + theme(legend.position = 'none')
  p
}