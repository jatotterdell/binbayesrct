#' Run Bayesian response-adaptive randomised trial
#'
#' @param n_seq        Sequence of interim analyses
#' @param dat          Potential outcomes data from `gen_potential_outcomes`
#' @param rand_quant   Quantity to decide allocations, `best` or `superior`
#' @param active_quant Quantity to decide active arms, `best` or `superior`
#' @param control      How to allocate to control, `fix`, `constant`, `match`.
#' @param approx       Use approximate rather than numerical methods.
#' @param perm_drop    Permanently drop any arm which is not `active`.
#' @param mwu          Use mass-weighted urn randomisation
#' @param minimum      Target minimum instead of maximum outcome
#' @param bpar         Parameter for control = `match`.
#' @param hpar1        Parameter to control amount of RAR
#' @param hpar2        Parameter to control scaling of RAR
#' @param h            Function for RAR
#' @param fpar1        Parameter to control lack of effect boundary (futility)
#' @param fpar2        Parameter to control lack of effect boundary (futility)
#' @param f            Function for utility boundary
#' @param ...          Other function arguments
#' @importFrom stats aggregate
#' @importFrom loo nlist
#' @export
beta_brar_trial <- function(
  n_seq,
  dat,
  rand_quant = "best",
  active_quant = "best",
  control = "fix",
  approx = FALSE,
  perm_drop = FALSE,
  mwu = FALSE,
  minimum = TRUE,
  bpar = 0.5,
  hpar1 = 2.5,
  hpar2 = 3,
  h = function(i) hpar1 * (i / length(n_seq)) ^ hpar2,
  fpar1 = 0.1,
  fpar2 = 1.5,
  f = function(i) fpar1 * (i / length(n_seq)) ^ fpar2,
  ...
) {
  if(max(n_seq) != nrow(dat)) stop("dimension mismatch, max(n_seq) should equal nrow(dat)")
  tdat <- dat[, 1:2]
  dat <- dat[, -(1:2)]

  n_max <- max(n_seq)
  n_int <- length(n_seq)
  n_seq_aug <- c(0, n_seq)
  n_new <- diff(n_seq_aug)

  A <- ncol(dat)
  a <- b <- 1 # Prior parameters

  n <- matrix(0, n_int + 1, A, dimnames = list("interim" = 0:n_int, "arm" = 1:A))
  y <- matrix(0, n_int + 1, A, dimnames = list("interim" = 0:n_int, "arm" = 1:A))
  p_supr <- matrix(0, n_int + 1, A - 1, dimnames = list("interim" = 0:n_int, "arm" = 2:A))
  p_best <- matrix(0, n_int + 1, A, dimnames = list("interim" = 0:n_int, "arm" = 1:A))
  p_best_trt <- matrix(0, n_int + 1, A - 1, dimnames = list("interim" = 0:n_int, "arm" = 2:A))
  p_best_active <- matrix(0, n_int + 1, A - 1, dimnames = list("interim" = 0:n_int, "arm" = 2:A))
  p_rand <- matrix(0, n_int + 1, A, dimnames = list("interim" = 0:n_int, "arm" = 1:A))

  is_active <- matrix(TRUE, n_int + 1, A, dimnames = list("interim" = 0:n_int, "arm" = 1:A))

  arm_mean <- matrix(0, n_int + 1, A, dimnames = list("interim" = 0:n_int, "arm" = 1:A))
  arm_var <- matrix(0, n_int + 1, A, dimnames = list("interim" = 0:n_int, "arm" = 1:A))
  eff_mean <- matrix(0, n_int + 1, A - 1, dimnames = list("interim" = 0:n_int, "arm" = 2:A))
  eff_var <- matrix(0, n_int + 1, A - 1, dimnames = list("interim" = 0:n_int, "arm" = 2:A))
  x <- factor(numeric(n_max), levels = 1:A)

  # Loop interim analyses
  for(i in 1:(n_int + 1)) {
    arm_mean[i, ] <- beta_mean(a + y[i, ], b + n[i, ] - y[i, ])
    arm_var[i, ] <- beta_var(a + y[i, ], b + n[i, ] - y[i, ])
    eff_mean[i, ] <- diff_beta_mean(a + y[i, ], b + n[i, ] - y[i, ])
    eff_var[i, ] <- diff_beta_var(a + y[i, ], b + n[i, ] - y[i, ])
    p_supr[i, ] <- beta_prob_supr(a + y[i, ], b + n[i, ] - y[i, ], reverse = minimum, approx = approx)
    p_best[i, ] <- beta_prob_best(a + y[i, ], b + n[i, ] - y[i, ], minimum = minimum, approx = approx)

    # Randomise on superiority or best?
    if(rand_quant == "superior") {
      q_rand = c(0, p_supr[i, ])
    } else if (rand_quant == "best") {
      q_rand = p_best[i, ]
    }

    # Drop on superiority or best?
    if(active_quant == "superior") {
      q_act = c(0, p_supr[i, ])
    } else if (active_quant == "best") {
      q_act = p_best[i, ]
    }

    # Active arm
    is_active[i, -1] <- q_act[-1] >= f(i - 1)
    if(perm_drop) {
      if(i > 1) {
        is_active[i, ] <- is_active[i, ] & is_active[i - 1, ]
      }
    }
    p_best_active[i, is_active[i, -1]] <- beta_prob_best(
      a + y[i, -1][is_active[i, -1]],
      b + n[i, -1][is_active[i, -1]] - y[i, -1][is_active[i, -1]],
      minimum = TRUE, approx = approx)

    # Fix control or match to best?
    if(control == "fix") {
      p_rand[i, ] <- fix_ctrl_brar(q_rand, is_active[i, ], h(i - 1))
    } else if(control == "constant") {
      p_rand[i, ] <- const_ctrl_brar(q_rand, is_active[i, ], h(i - 1))
    } else {
      p_rand[i, ] <- match_ctrl_brar(q_rand, is_active[i, ], n[i, ], h(i - 1), bpar)
    }

    if(i < n_int + 1) {
      if(mwu) {
        new_x <- factor(mass_weighted_urn_design(p_rand[i, ], n_new[i], alpha = 5)$trt, levels = 1:A)
      } else {
        new_x <- factor(sample.int(A, n_new[i], prob = p_rand[i, ], replace = T), levels = 1:A)
      }

      new_n <- table(new_x, dnn = NULL)
      current_indx <- (n_seq_aug[i] + 1):n_seq[i]
      new_y <- aggregate(dat[cbind(current_indx, new_x)], list(new_x), sum, drop = F)[, 2]
      new_y[is.na(new_y)] <- 0

      x[current_indx] <- new_x
      n[i + 1, ] <- n[i, ] + new_n
      y[i + 1, ] <- y[i, ] + new_y
    }
  }
  return(nlist(
    n, y,
    p_best, p_supr, p_rand, p_best_active,
    arm_mean, arm_var, eff_mean, eff_var, is_active
  ))
}
