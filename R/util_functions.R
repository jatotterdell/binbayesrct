#' Basic Moore-Penrose pseudo-inverse
#'
#' @param X A matrix
#' @return Pseudo-inverse of X
#' @export
psolve <- function(X) {
  solve(crossprod(X), t(X))
}


#' Normalise vector to sum to 1
#'
#' @param x A numeric vector
#' @return The vector `x` normalised to sum to 1.
normalise <- function(x)
  x / sum(x)


#' Find first element satisfying condition
#'
#' @param x A logical vector
#' @param v A default value if none satisfy condition
#' @return Index of first to satisfy condition
#' @export
findfirst <- function(x, v = NA) {
  j <- which(x)
  if(length(j)) min(j) else v
}

#' Compare Normal approximation to Beta
#'
#' @param a First Beta parameter
#' @param b Second Beta parameter
#' @param ... Other arguments to `curve()`
#' @return A plot of Beta density and Normal approximation
#'
#' @export
plot_beta_norm <- function(a, b, ...) {
  fbeta <- function(x) dbeta(x, a, b)
  fnorm <- function(x) dnorm(x, a/(a + b), sqrt( a*b / ((a+b)^2*(a+b+1))))
  graphics::curve(fbeta, ...)
  graphics::curve(fnorm, add = TRUE, col = "red", ...)
}

#' Calculate density of beta-binomial distribution
#'
#' @param x The value at which to evaluate the density
#' @param n The sample size
#' @param a First parameter
#' @param b Second parameter
#'
#' @return Value of beta-binomial(n,a,b) evaluated at x
#'
#' @examples
#' dbetabinom(5, 10, 2, 3)
#'
#' @export
dbetabinom <- function(x, n, a = 1, b = 1){
  if(!(all(c(a, b) > 0))) stop("a and b must be > 0")
  if(any(n < 1)) stop("n must be > 0")
  if(any(x < 0)) stop("x must be >= 0")

  num <- lgamma(a + b) + lgamma(n + 1) + lgamma(x + a) + lgamma(n - x + b)
  den <- lgamma(a) + lgamma(b) + lgamma(x + 1) + lgamma(n - x + 1) + lgamma(n + a + b)
  prob <- exp(num - den)
  prob
}

#' Draw random variates from beta-binomial distribution
#'
#' @import stats
#'
#' @param n The number of random values to sample
#' @param m The sample size
#' @param a First parameter
#' @param b Second parameter
#'
#' @examples
#' rbetabinom(2, 10, 2, 3)
#'
#' @export
rbetabinom <- function(n, m, a = 1, b = 1) {
  if(!(all(c(a, b) > 0))) stop("a and b must be > 0")
  if(!(all(n > 0))) stop("n must be > 0")

  stats::rbinom(n, m, stats::rbeta(n, a, b))
}

#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using numerical integration
#'
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#' @param ... other arguments passed to integrate/quadgk function
#'
#' @return The value of the integral
#'
#' @examples
#' beta_ineq(5, 5, 3, 7)
#'
#' @export
beta_ineq <- function(a, b, c, d, delta = 0, ...) {

  if(!(all(c(a, b, c, d) > 0))) stop("a, b, c, d must be > 0")

  integrand <- function(x) { stats::dbeta(x, a, b)*stats::pbeta(x - delta, c, d) }
  tryCatch(
    pracma::quadgk(integrand, delta, 1, ...),
    error = function(err) NA)
}

#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using Normal approximation.
#'
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#' @return The value of the integral
#'
#' @examples
#' beta_ineq_approx(5, 5, 3, 7)
#'
#' @export
beta_ineq_approx <- function(a, b, c, d, delta = 0) {
  if(!(all(c(a, b, c, d) > 0))) stop("a, b, c, d must be > 0")

  m1 <- a / (a + b)
  v1 <- a*b / ( (a + b)^2 * (a + b + 1))
  m2 <- c / (c + d)
  v2 <- c*d / ( (c + d)^2 * (c + d + 1))
  z <- (m1 - m2 - delta) / sqrt(v1 + v2)
  return(stats::pnorm(z))
}

#' Calculate Pr(X > Y + delta) where X and Y are independent Beta random variables
#' using Monte Carlo method.
#'
#' @param a Parameter one of beta density for X
#' @param b Parameter two of beta density for X
#' @param c Parameter one of beta density for Y
#' @param d Parameter two of beta density for Y
#' @param delta The difference we wish to assess (i.e. X - Y > delta)
#' @param sims The number of Monte Carlo variates to generate for estimation
#'
#' @return The value of the integral
#'
#' @examples
#' beta_ineq_sim(5, 5, 3, 7)
#'
#' @export
beta_ineq_sim <- function(a, b, c, d, delta = 0, sims = 10000) {
  if(!(all(c(a, b, c, d) > 0))) stop("a, b, c, d must be > 0")

  lens <- unlist(lapply(list(a, b, c, d), length))
  if(any(max(lens)- min(lens) != 0)) stop("a, b, c, d must be same len")

  X <- lapply(1:length(a), function(x) stats::rbeta(sims, a[x], b[x]))
  Y <- lapply(1:length(a), function(x) stats::rbeta(sims, c[x], d[x]))

  means <- lapply(1:length(a), function(x) mean(X[[x]] > Y[[x]] + delta))
  unlist(means)
}



N <- NULL
P <- NULL
#' Calculate the predicted probability of success
#'
#' @import data.table
#'
#' @param a First parameter of first beta random variable
#' @param b Second parameter of first beta random variable
#' @param c First paramter of second beta random variable
#' @param d Second parameter of second beta random variable
#' @param m1 Sample size to predict for first beta random variable
#' @param m2 Sample size to predict for second beta random variable
#' @param k_ppos The posterior probability cut-point to be assessed
#' @param post_method Method for calcuation: exact, approx, or sim.
#' @param post_sim Number of posterior draws to use for PPoS
#'
#' @return The predicted probability of success
#'
#' @export
calc_ppos <- function(
  a,
  b,
  c,
  d,
  m1,
  m2,
  k_ppos,
  post_method = "exact",
  post_sim = 1e4) {

  requireNamespace(data.table)
  if(!(all(c(a, b, c, d, m1, m2) > 0))) stop("a, b, c, d, m1, m2 must be > 0")
  if(k_ppos < 0 | k_ppos > 1) stop("k_ppos must be in [0, 1]")

  calc_post <- switch(post_method,
                      "exact" = beta_ineq,
                      "approx" = beta_ineq_approx,
                      "sim" = beta_ineq_sim)

  y1pred <- rbetabinom(post_sim, m1, a, b)
  y2pred <- rbetabinom(post_sim, m2, c, d)
  ypred <- data.table::data.table(y1pred = y1pred, y2pred = y2pred)[, .N, keyby = list(y1pred, y2pred)]
  ypred[, data.table::`:=`(P = Vectorize(calc_post)(a + y1pred,
                                        b + m1 - y1pred,
                                        c + y2pred,
                                        d + m2 - y2pred))]
  ypred[, c(sum(N*(P > k_ppos)) / sum(N))]
}


#' Mass-weighted urn randomisation
#'
#' @param target_alloc The target allocation ratios
#' @param sample_size The number of allocations to generate
#' @param alpha Parameter to control imbalance between arms
#' @return A list detailing the mass-weighted-urn process.
#' @export
#' @importFrom stats runif
mass_weighted_urn_design <- function(
  target_alloc,
  sample_size,
  alpha = 4
) {
  arms <- length(target_alloc)
  prob_alloc <- target_alloc / sum(target_alloc)
  # Masses
  x <- matrix(0, sample_size + 1, arms)
  x[1, ] <- alpha * prob_alloc
  # Sample size
  n <- matrix(0, sample_size + 1, arms)
  # Random number
  y <- runif(sample_size)
  # Conditional selection probability
  p <- matrix(0, sample_size + 1, arms)
  # Imbalance
  d <- rep(0, sample_size)
  # Allocation Predictability
  g <- rep(0, sample_size + 1)
  # Treatment assignment
  trt <- rep(0, sample_size)

  imbalance_cap <- sqrt(sum(((alpha - 1)*(1 - prob_alloc) + (arms - 1))^2))

  for(i in 2:(sample_size + 1)) {
    # Update allocation probabilities
    p[i - 1, ] <- pmax(alpha * prob_alloc - n[i - 1, ] + (i - 1)*prob_alloc, 0)
    p[i - 1, ] <- p[i - 1, ] / sum(p[i - 1, ])
    trt[i-1] <- findInterval(y[i - 1], c(0, cumsum(p[i - 1, ])))
    # Update sample sizes
    n[i, ] <- n[i - 1, ]
    n[i, trt[i-1]] <- n[i, trt[i-1]] + 1
    # Update urn masses
    x[i, trt[i-1]] <- x[i - 1, trt[i-1]] - 1 + prob_alloc[trt[i-1]]
    x[i, -trt[i-1]] <- x[i - 1, -trt[i-1]] + prob_alloc[-trt[i-1]]
    # Calculate imbalance
    d[i - 1] <- sqrt(sum((n[i, ] - (i - 1)*prob_alloc)^2))
    # Calculate allocation predictability
    g[i] <- d[i - 1] / alpha
  }
  return(list(
    max_imbalance_bound = imbalance_cap,
    imbalance = d,
    alloc_predict = g,
    rand_num = y,
    trt = trt,
    mass = x,
    sample_size = n,
    selection_prob = p))
}


#' Mean of a beta random variable with parameters a and b
#'
#' @param a Par 1
#' @param b Par 2
#' @examples
#' beta_mean(5, 5)
#' @export
beta_mean <- function(a, b) {
  return(a / (a + b))
}


#' Mean of added beta random variables
#'
#' Always taken in reference to control, e.g. `(a[1], b[1])`.
#'
#' @param a Collection of shape 1 par
#' @param b Collection of shape 2 par
#' @examples
#' diff_beta_mean(c(5, 5, 4), c(5, 5, 4))
#' @export
diff_beta_mean <- function(a, b) {
  return(beta_mean(a[-1], b[-1]) - beta_mean(a[1], b[1]))
}


#' Variance of a beta random variable with parameters a and b
#'
#' @param a Par 1
#' @param b Par 2
#' @examples
#' beta_var(5, 5)
#' @export
beta_var <- function(a, b) {
  return( a * b / ( (a + b)^2 * (a + b + 1)) )
}


#' Variance of a added beta random variables
#'
#' Always taken in reference to control, e.g. `(a[1], b[1])`.
#'
#' @param a Collection of shape 1 par
#' @param b Collection of shape 2 par
#' @examples
#' diff_beta_var(c(5, 5, 4), c(5, 5, 4))
#' @export
diff_beta_var <- function(a, b) {
  return(beta_var(a[1], b[1]) + beta_var(a[-1], b[-1]))
}


#' Probability that each arm superior to reference.
#'
#' Calculates event that one Beta RV larger (or smaller if reverse) to another Beta RV by integration.
#'
#' @param a First value is reference
#' @param b First value is reference
#' @param reverse Reverse direction of superiority
#' @param ... Other arguments to `pracma::quadgk`
#' @examples
#' beta_prob_supr_numeric(c(5, 10, 15), c(15, 10, 5))
#' beta_prob_supr_numeric(c(5, 10, 15), c(15, 10, 5), reverse = TRUE)
#' @importFrom pracma quadgk
#' @importFrom stats dbeta pbeta
#' @export
beta_prob_supr_numeric <- function(a, b, reverse = FALSE, ...) {
  if(!(all(c(a, b) > 0))) stop("a, b, must be > 0")
  k <- length(a)
  ans <- numeric(k-1)
  for(i in 2:k) {
    f <- function(z)  dbeta(z, a[1], b[1])*(1 - pbeta(z, a[i], b[i]))
    ans[i-1] <- tryCatch(
      pracma::quadgk(f, 0, 1, ...),
      error = function(err) NA)
  }
  if(reverse) ans <- 1 - ans
  return(ans)
}


#' Probability that each arm superior to reference.
#'
#' Calculates event that one Beta RV larger (or smaller if reverse) to another Beta RV by normal approximation.
#'
#' @param a First value is reference
#' @param b First value is reference
#' @param reverse Reverse direction of superiority
#' @examples
#' beta_prob_supr_approx(c(5, 10, 15), c(15, 10, 5))
#' beta_prob_supr_approx(c(5, 10, 15), c(15, 10, 5), reverse = TRUE)
#' @importFrom stats pnorm
#' @export
beta_prob_supr_approx <- function(a, b, reverse = FALSE) {
  if(!(all(c(a, b) > 0))) stop("a, b, must be > 0")
  m <- a / (a + b)
  v <- a*b / ( (a + b)^2 * (a + b + 1))
  z <- (m[-1] - m[1]) / sqrt(v[-1] + v[1])
  ans <- pnorm(z)
  if(reverse) ans <- 1 - ans
  return(ans)
}


#' Probability that each arm superior to reference.
#'
#' Calculates event that one Beta RV larger (or smaller) to another Beta RV by normal approximation.
#'
#' @param a First value is reference
#' @param b First value is reference
#' @param reverse Reverse direction of superiority
#' @param approx Use normal approximation instead of numerical integration
#' @param ... Other arguments to `pracma::quadgk`
#' @examples
#' beta_prob_supr(c(5, 10, 15), c(15, 10, 5))
#' beta_prob_supr(c(5, 10, 15), c(15, 10, 5), approx = TRUE)
#' @export
beta_prob_supr <- function(a, b, reverse = FALSE, approx = FALSE, ...) {
  if(!(all(c(a, b) > 0))) stop("a, b, must be > 0")
  k <- length(a)
  ans <- numeric(k-1)
  if(!approx) {
    return(beta_prob_supr_numeric(a, b, reverse, ...))
  } else {
    return(beta_prob_supr_approx(a, b, reverse))
  }
}


#' Probability that each element of a collection of Beta RV is the maximum/minimum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum Use minimum instead of maximum
#' @param ... Other arguments to `pracma::quadgk`
#' @return A vector of `length(a)` giving `Pr(a[i] is max)`
#' @examples
#' beta_prob_best_numeric(1 + seq(1, 9, 2), 1 + rep(10, 5))
#' @export
beta_prob_best_numeric <- function(a, b, minimum = FALSE, ...) {
  k <- length(a)
  ans <- numeric(k)
  for(i in 1:k) {
    if(!minimum) {
      f <- function(z)
        sapply(z, function(x) dbeta(x, a[i], b[i]) * prod(pbeta(x, a[-i], b[-i])))
    } else {
      f <- function(z)
        sapply(z, function(x) dbeta(x, a[i], b[i]) * prod(1 - pbeta(x, a[-i], b[-i])))
    }
    ans[i] =   tryCatch(
      pracma::quadgk(f, 0, 1, ...),
      error = function(err) NA)
  }
  return(ans)
}


#' Probability that each element of a collection of Beta RV is the maximum/minimum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum Use minimum instead of maximum
#' @param nsim Number of `rbeta` draws to use in assessment.
#' @return A vector of `length(a)` giving `Pr(a[i] is max)`
#' @examples
#' beta_prob_best_approx(1 + seq(1, 9, 2), 1 + rep(10, 5))
#' @importFrom stats rbeta
#' @export
beta_prob_best_approx <- function(a, b, minimum = FALSE, nsim = 2e4) {
  k <- length(a)
  ans <- numeric(k)
  r <- matrix(rbeta(k*nsim, a, b), nsim, k, byrow = T)
  if(!minimum) {
    ans <- prop.table(table(factor(max.col(r), levels = 1:k, labels = 2:(k+1))))
  } else {
    ans <- prop.table(table(factor(max.col(1 - r), levels = 1:k, labels = 2:(k+1))))
  }
  return(ans)
}


#' Probability that each element of a collection of Beta RV is the maximum/minimum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum Use minimum instead of maximum
#' @param approx Use approximation instead of integration
#' @param ... Other arguments to `pracma::quadgk`
#' @return A vector of `length(a)` giving `Pr(a[i] is max)`
#' @export
beta_prob_best <- function(a, b, minimum = FALSE, approx = FALSE, ...) {
  if(!(all(c(a, b) > 0))) stop("a, b, must be > 0")
  k <- length(a)
  ans <- numeric(k)
  if(!approx) {
    return(beta_prob_best_numeric(a, b, minimum, ...))
  } else {
    return(beta_prob_best_approx(a, b, minimum))
  }
}


#' Entropy of the response density of the unknown maximum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum If true use minimum instead of maximum.
#' @param ... Other arguments to `quadgk`.
#' @importFrom pracma quadgk
ent_best_arm_numeric <- function(a, b, minimum = FALSE, ...) {
  k <- length(a)
  if(!minimum) {
    f <- function(x) {
      ss <- sum(sapply(1:k, function(i) {
        dbeta(x, a[i], b[i]) * prod(pbeta(x, a[-i], b[-i]))
      }))
      ss <- ss * log(ss)
      if(is.na(ss))
        ss <- 0
      return(-ss)
    }
  } else {
    f <- function(x) {
      ss <- sum(sapply(1:k, function(i) {
        dbeta(x, a[i], b[i]) * prod(1 - pbeta(x, a[-i], b[-i]))
      }))
      ss <- ss * log(ss)
      if(is.na(ss))
        ss <- 0
      return(-ss)
    }
  }
  ans <- pracma::quadgk(Vectorize(f), 0, 1, ...)
  return(ans)
}


#' Evaluate joint density of maximum of Beta RVs
#'
#' @param x Value to evaluate density at
#' @param a Vector of shape1 parameters
#' @param b Vector of shape1=2 parameters
#' @param minimum Use minimum instead of maximum.
#' @export
beta_best_dens <- function(x, a, b, minimum = FALSE) {
  k <- length(a)
  if(!minimum) {
    sapply(x, function(z) sum(sapply(1:k, function(i) dbeta(z, a[i], b[i]) * prod(pbeta(z, a[-i], b[-i])))))
  } else {
    sapply(x, function(z) sum(sapply(1:k, function(i) dbeta(z, a[i], b[i]) * prod(1 - pbeta(z, a[-i], b[-i])))))
  }
}


#' Entropy of the response density of the unknown maximum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum If true use minimum instead of maximum.
#' @param nn Grid points.
#' @importFrom pracma quadgk
#' @export
ent_best_arm_approx <- function(a, b, minimum = FALSE, nn = 100) {
  k <- length(a)
  del <- 1/nn
  x <- seq(0 + del, 1 - del, length.out = nn - 1)
  f <- beta_best_dens(x, a, b, minimum)
  f[f == 0] <- .Machine$double.eps^12
  return(-sum(f*log(f) * del))
}


#' Entropy of the response density of the unknown maximum
#'
#' @param a Vector of shape1 parameters
#' @param b Vector of shape2 parameters
#' @param minimum If true use minimum instead of maximum.
#' @param approx Use approximation instead of numerical.
#' @param ... Other arguments to `quadgk`.
#' @export
ent_best_arm <- function(a, b, minimum = F, approx = F, ...) {
  if(!(all(c(a, b) > 0))) stop("a, b, must be > 0")
  if(!approx) {
    return(ent_best_arm_numeric(a, b, minimum, ...))
  } else {
    return(ent_best_arm_approx(a, b, minimum))
  }
}


#' Calculate change in utility using variance
#'
#' @param n Sample size in each arm
#' @param y Response in each arm
#' @param a Prior shape
#' @param b Prior shape
#' @export
utility_var <- function(n, y, a = 1, b = 1) {
  K <- length(n)
  # Prior variance
  priorV <- diff_beta_var(a, b)
  priorU <- sum(priorV)
  # Current utility
  V <- diff_beta_var(y + a, n - y + b)
  U <- sum(V)
  # Possible allocations
  new_n <- t(n + diag(1, K))
  # If response to observed from allocation
  new_y <- t(y + diag(1, K))
  # For each arm, what is the expected gain in utility
  Eu <- numeric(K)
  for(j in 1:K) {
    # predicted probability of event
    m <- (a + y[j]) / (a + b + n[j])
    # utility if event occurs
    H1 <- diff_beta_var(new_y[j, ] + a, new_n[j, ] - new_y[j, ] + b)
    # utility if no event occurs
    H0 <- diff_beta_var(y + a, new_n[j, ] - y + b)
    # Expected utility
    Eu[j] <- m * sum(H1) + (1 - m) * sum(H0)
  }
  # Change in utility
  return(U - Eu)
}


#' Calculate change in utility using entropy
#'
#' @param n Sample size in each arm
#' @param y Response in each arm
#' @param a Prior shape
#' @param b Prior shape
#' @param ... Other function arguments
#' @examples
#' utility_ent(rep(0, 4), rep(0, 4))
#' utility_ent(c(20, 19, 18, 17), c(8,9,10,11))
#' @export
utility_ent <- function(n, y, a = 1, b = 1, ...) {
  K <- length(n)
  # Current utility
  H <- ent_best_arm(y + a, n - y + b, ...)
  # Possible allocations
  new_n <- t(n + diag(1, K))
  # If response to observed from allocation
  new_y <- t(y + diag(1, K))
  # For each arm, what is the expected gain in utility
  Eu <- numeric(K)
  for(j in 1:K) {
    # predicted probability of event
    m <- (a + y[j]) / (a + b + n[j])
    # utility if event occurs
    H1 <- ent_best_arm(new_y[j, ] + a, new_n[j, ] - new_y[j, ] + b, ...)
    # utility if no event occurs
    H0 <- ent_best_arm(y + a, new_n[j, ] - y + b, ...)
    # Expected utility
    Eu[j] <- m * H1 + (1 - m) * H0
  }
  # Change in utility
  return(H - Eu)
}



#' Generate potential binary outcomes
#'
#' @param n Sample size
#' @param p Vector of probabilities of outcome
#' @param accrual_rate A function target integer argument which
#' generates that number of arrival times.
#' @param delay The time until endpoint is observed.
#' @param ... Optional additional arguments to the accrual_rate function.
#' @return A matrix of potential outcomes for each arm
#' @importFrom stats rbinom
#' @export
gen_potential_outcomes <- function(
  n,
  p,
  accrual_rate = function(n) (1:n)/2,
  delay = 0,
  ...) {

  k <- length(p)
  tim <- sort(accrual_rate(n, ...))
  obs <- tim + delay
  cbind(
    t_acc = tim,
    t_obs = obs,
    matrix(
      rbinom(n*k, 1, p), n, k, byrow = T,
      dimnames = list('id' = 1:n, 'o' = paste0('y_', 1:k))
    ))
}


#' Probability each column is maximum (or minimum)
#'
#' @param mat A matrix of draws.
#' @param minimum If `TRUE` determine best based on minimum instead of maximum value.
#' @export
prob_best <- function(mat, minimum = F) {
  if(!minimum) {
    as.numeric(prop.table(table(factor(max.col(mat), levels = 1:ncol(mat)))))
  } else {
    as.numeric(prop.table(table(factor(max.col(-mat), levels = 1:ncol(mat)))))
  }
}


#' BRAR with fixed ratio to control
#'
#' @param quantity The quantity for allocation ratios
#' @param active Flag for active arms
#' @param h Scaling factor
#' @return A vector of allocation probabilities
#' @export
fix_ctrl_brar <- function(quantity, active, h) {
  arms <- length(active)
  w <- numeric(arms)
  if(!any(active)) {
    w <- rep(0, arms)
  } else if(!any(active[-1])) {
    w[1] <- 1
  } else if(!active[1]) {
    w[1] <- 0
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
  } else {
    w[1] <- 1 / (sum(active))
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
    w[-1] <- (1 - w[1]) * w[-1]
  }
  return(w)
}


#' BRAR with fixed ratio to control
#'
#' @param quantity The quantity for allocation ratios
#' @param active Flag for active arms
#' @param h Scaling factor
#' @param fixq The quantity to fix control at
#' @export
fix_ctrl_brar2 <- function(quantity, active, h, fixq) {
  arms <- length(active)
  w <- numeric(arms)
  if(!any(active)) {
    w <- rep(0, arms)
  } else if(!any(active[-1])) {
    w[1] <- 1
  } else if(!active[1]) {
    w[1] <- 0
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
  } else {
    w[1] <- fixq
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
    w[-1] <- (1 - w[1]) * w[-1]
  }
  return(w)
}


#' BRAR with matching control on biggest sample size (Trippa)
#'
#' @param quantity The quantity for allocation ratios
#' @param active Flag for active arms
#' @param n Sample size in each arm
#' @param h Scaling factor
#' @param bpar Control matching parameter
match_ctrl_brar <- function(quantity, active, n, h, bpar = 1) {
  arms <- length(active) + 1
  w <- numeric(arms)
  if(!any(active)) {
    w[1] <- 1
    return(w)
  } else {
    w[-1] <- quantity ^ h
    w[-1][!active] <- 0
    C <- sum(w[-1]) / sum(active)
    w[1] <- C * exp(bpar * (max(n[-1]) - n[1]))
    w <- normalise(w)
    return(w)
  }
}


#' BRAR with constant proportion to control
#'
#' @param quantity The quantity to be used for allocation ratios.
#' @param active Flag indicating active arms.
#' @param h Scaling factor value
#' @return Allocation probabilities to each arm.
#' @examples
#' const_ctrl_brar(runif(5), sample(c(TRUE, FALSE), 5, rep = TRUE), 1)
#' @export
const_ctrl_brar <- function(quantity, active, h) {
  arms <- length(active)
  w <- numeric(arms)
  if(!any(active)) {
    w <- rep(0, arms)
  } else if(!active[1]) {
    w[1] <- 0
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
  } else {
    w[1] <- 1 / arms
    w[-1] <- quantity[-1] ^ h
    w[-1][!active[-1]] <- 0
    w[-1] <- normalise(w[-1])
    w[-1] <- (1 - w[1]) * w[-1]
  }
  return(w)
}

#' BRAR on all arms.
#'
#' @param quantity The quantity to be used for allocation ratios.
#' @param active Flag indicating active arms.
#' @param h Scaling factor value
#' @return Allocation probabilities to each arm.
#' @examples
#' brar_all(runif(5), sample(c(TRUE, FALSE), 5, rep = TRUE), 1)
#' @export
brar_all <- function(quantity, active, h) {
  arms <- length(active)
  w <- numeric(arms)
  if(!any(active)) {
    w <- rep(0, arms)
  } else {
    w <- quantity ^ h
    w[!active] <- 0
    w <- normalise(w)
  }
  return(w)
}
