# =========== CI dispatcher (YOU fill these) ===========

ci_dispatch <- list(
  cox        = function(obj) stop("CI fn for 'cox' not implemented"),
  flexsurv = function(obj) stop("CI fn for 'flexsurv' not implemented"),
  msm        = function(obj) stop("CI fn for 'msm' not implemented"),
  msm_age    = function(obj) stop("CI fn for 'msm_age' not implemented"),
  nhm        = function(obj) stop("CI fn for 'nhm' not implemented"),
  mipd_iter  = function(obj) stop("CI fn for 'mipd_iter' not implemented")

)
# ======================================================



# ----------------- helpers -----------------
parse_seed <- function(fname) {
  b <- basename(fname)
  m <- str_match(b, "(?i)seed[_-]?([0-9]+)")
  if (!is.na(m[1,2])) return(as.integer(m[1,2]))
  NA_integer_
}

safe_readRDS <- purrr::safely(readRDS, otherwise = NULL)

safe_ci_call <- function(tag, obj) {
  fn <- ci_dispatch[[tag]]
  if (!is.function(fn))
    stop(sprintf("No CI function provided for tag '%s'", tag))
  out <- fn(obj)
  req <- c("param","est","lwr","upr")
  if (!all(req %in% names(out)))
    stop(sprintf("CI fn for '%s' must return columns: %s",
                 tag, paste(req, collapse=", ")))
  tibble::as_tibble(out) %>% mutate(param = as.character(.data$param))
}

summarize_across_seeds <- function(per_seed) {
  per_seed %>%
    group_by(scheme_id, model, param, trans_idx) %>%
    summarise(
      n_seeds         = n_distinct(seed),
      est_mean        = mean(est, na.rm = TRUE),
      est_sd          = sd(est, na.rm = TRUE),
      rel_bias        = mean(rel_bias, na.rm =TRUE),
      coverage_rate   = mean(coverage, na.rm = TRUE),
      type_I_rate      = mean(typeI, na.rm = TRUE),
      type_II_rate     = mean(typeII, na.rm = TRUE)
    ) %>%
    arrange(param, trans_idx, model, scheme_id)
}



list_model_files <- function(sim_root, scheme_id, term_rate, tag) {
  base_dir <- file.path(sim_root, sprintf("cov_scheme_%d", scheme_id), term_rate, tag)

  # helper to list files by regex + attach seed + submodel label
  .collect <- function(dir, pattern, sub = NA_character_) {
    if (!dir.exists(dir)) return(tibble::tibble(filepath=character(), seed=integer(), model_sub=character()))
    fs_all <- list.files(dir, pattern="\\.[Rr][Dd][Ss]$", full.names=TRUE, recursive=FALSE)
    fs <- fs_all[grepl(pattern, basename(fs_all))]
    if (!length(fs)) return(tibble::tibble(filepath=character(), seed=integer(), model_sub=character()))
    tibble::tibble(
      filepath = fs,
      seed     = vapply(fs, parse_seed, integer(1), USE.NAMES = FALSE),
      model_sub= sub
    ) |>
      dplyr::mutate(seed = dplyr::if_else(is.na(seed), dplyr::row_number(), seed)) |>
      dplyr::arrange(seed)
  }

  if (tag == "cvae1") {
    dir <- file.path(base_dir, "fits")
    return(.collect(dir, pattern = ".", sub = tag)) # keep sub as tag for clarity
  } else if (tag == "flexsurv") {
    dir <- base_dir
    df_panel <- .collect(dir, pattern = "^seed_[0-9]+_model_flex_panel\\.rds$", sub = "flexsurv_panel")
    df_real  <- .collect(dir, pattern = "^seed_[0-9]+_model_flex_real\\.rds$",  sub = "flexsurv_real")
    return(dplyr::bind_rows(df_panel, df_real))
  } else if (tag == "mipd_iter") {
    dir <- base_dir
    pat_model <- "^seed_[0-9]+_model_mipd\\.[Rr][Dd][Ss]$"
    pat_bias  <- "^seed_[0-9]+_bias_mipd\\.[Rr][Dd][Ss]$"
    df_model <- .collect(dir, pattern = pat_model, sub = "mipd_iter_model") %>%
      dplyr::mutate(kind = "model")
    df_bias  <- .collect(dir, pattern = pat_bias,  sub = "mipd_iter_bias") %>%
      dplyr::mutate(kind = "bias")
    return(dplyr::bind_rows(df_model, df_bias))
  }else {
    dir <- base_dir
    pat <- paste0("^seed_[0-9]+_model_", tag, "\\.[Rr][Dd][Ss]$")
    return(.collect(dir, pattern = pat, sub = tag))
  }
}
# =============== ci per model type ============================

# coverage = 1 when the true value is inside estimated ic
# then for covs when feasible
# type I error = 1 when true value is not significant but estimated as significant (feasible when true value zero)
# type II error = 1  when true value significant but estimated as not significant (feasible when true value !zero)

ci_dispatch$cox <- function(obj) {

  ci  <- suppressMessages(suppressWarnings(stats::confint(obj)))
  est <- stats::coef(obj)


  df <- tibble::tibble(
    trans_idx = seq_along(est),
    param     = "cov",
    est       = as.numeric(est),
    lwr       = ci[, 1],
    upr       = ci[, 2]
  )


  gt <- as.matrix(ground_truth)
  j  <- match(df$param, colnames(gt))
  ok <- !is.na(j) & df$trans_idx >= 1 & df$trans_idx <= nrow(gt)

  gt_vals <- rep(NA_real_, nrow(df))
  gt_vals[ok] <- gt[cbind(df$trans_idx[ok], j[ok])]

  df$gt <- gt_vals
  df$rel_bias <- if_else(gt_vals != 0, (df$est- gt_vals) / abs(gt_vals), NA_real_)
  df$coverage <- as.integer( gt_vals >= df$lwr & gt_vals <= df$upr)

  df$typeI    <- ifelse( gt_vals == 0,
                         as.integer(df$lwr > 0 | df$upr < 0),
                         NA_integer_)

  df$typeII <- ifelse(abs(gt_vals) > 0,
                      as.integer(df$lwr < 0 & df$upr > 0), NA_integer_)

  df


}

ci_dispatch$flexsurv <- function(obj) {

  ci <- lapply(obj, function(model) confint(model))
  est <- lapply(obj, function(model) coef(model))


  df <- map_dfr(seq_along(ci), function(k) {
    ci_k   <- as.matrix(ci[[k]])
    pars   <- rownames(ci_k)
    est_k  <- est[[k]]
    est_k  <- setNames(as.numeric(est_k), names(est_k))

    tibble(
      trans_idx = k,
      param     = pars,
      est       = unname(est_k[pars]),
      lwr       = ci_k[, 1],
      upr       = ci_k[, 2]
    )
  })

  param_order <- c("shape", "rate", "cov")
  df <- df %>%
    mutate(param = factor(param, levels = param_order)) %>%
    arrange(param, trans_idx)

  gt <- as.matrix(ground_truth)
  j  <- match(df$param, colnames(gt))
  ok <- !is.na(j) & df$trans_idx >= 1 & df$trans_idx <= nrow(gt)

  gt_vals <- rep(NA_real_, nrow(df))
  gt_vals[ok] <- gt[cbind(df$trans_idx[ok], j[ok])]

  df$gt <- gt_vals
  df$rel_bias <- if_else(gt_vals != 0, (df$est- gt_vals) / abs(gt_vals), NA_real_)
  df$coverage <- as.integer(!is.na(gt_vals) & gt_vals >= df$lwr & gt_vals <= df$upr)

  df$typeI    <- ifelse( gt_vals == 0,
                         as.integer(df$lwr > 0 | df$upr < 0),
                         NA_integer_)

  df$typeII <- ifelse(abs(gt_vals) > 0,
                      as.integer(df$lwr < 0  &  df$upr > 0), NA_integer_)

  df


}

ci_dispatch$msm <- function(obj) {

  ci  <- obj$ci
  est <- obj$estimates

  df <- tibble::tibble(
    trans_idx = c(1,2,3),
    param     = "cov",
    est       = as.numeric(est[4:6]),
    lwr       = ci[4:6, 1],
    upr       = ci[4:6, 2]
  )

  gt <- as.matrix(ground_truth)
  j  <- match(df$param, colnames(gt))
  ok <- !is.na(j) & df$trans_idx >= 1 & df$trans_idx <= nrow(gt)

  gt_vals <- rep(NA_real_, nrow(df))
  gt_vals[ok] <- gt[cbind(df$trans_idx[ok], j[ok])]

  df$gt <- gt_vals
  df$rel_bias <- if_else(gt_vals != 0, (df$est- gt_vals) / abs(gt_vals), NA_real_)
  df$coverage <- as.integer( gt_vals >= df$lwr & gt_vals <= df$upr)

  df$typeI    <- ifelse( gt_vals == 0,
                         as.integer(df$lwr > 0 | df$upr < 0),
                         NA_integer_)

  df$typeII <- ifelse(abs(gt_vals) > 0,
                      as.integer(df$lwr < 0  &  df$upr > 0), NA_integer_)

  df


}
ci_dispatch$msm_age <- ci_dispatch$msm

# can be used for msm_age as well as long as we don't require estimates of the ic for
# rate and shape in that case we need bootstrap


get_params_nhm <- function (x, ci = TRUE, ...)
{
  model <- x
  par <- model$par
  fisher <- model$hess
  par_names <- model$parnames
  if (is.null(fisher))
    fisher <- model$fisher
  if (is.null(fisher) & ci)
    stop("Fisher information required for confidence intervals")
  if (is.null(par_names))
    par_names <- sapply(1:length(par), function(x) paste("Parameter",
                                                         x))
  if (model$singular & ci) {
    warning("Cannot calculate confidence intervals due to a singular Hessian")
    ci <- FALSE
  }
  if (ci) {
    std.err <- diag(solve(fisher))^0.5
    mat <- round(cbind(par, par - qnorm(0.975) * std.err,
                       par + qnorm(0.975) * std.err), 4)
    dimnames(mat)[[2]] <- c("Est", "Low 95%", "Up 95%")
  }
  else {
    mat <- round(cbind(par), 4)
    dimnames(mat)[[2]] <- c("Est")
  }
  if (!is.null(model$fixedpar)) {
    mat2 <- array(NA, c(model$npar, dim(mat)[2]))
    mat2[-model$fixedpar, ] <- mat
    mat2[model$fixedpar, 1] <- model$fixval
    dimnames(mat2)[[2]] <- dimnames(mat)[[2]]
    mat <- mat2
  }
  dimnames(mat)[[1]] <- par_names
  return(mat)

}

ci_dispatch$nhm <- function(obj) {

  ci  <- get_params_nhm(obj, ci = TRUE)
  ci <-  unname(ci[,2:3])
  est <- obj$par

  df <- tibble::tibble(
    trans_idx = rep(1:3,3),
    param     = c(rep("rate",3), rep("shape",3), rep("cov",3)),
    est       = as.numeric(est),
    lwr       = ci[, 1],
    upr       = ci[, 2]
  )


  param_order <- c("shape", "rate", "cov")
  df <- df %>%
    mutate(param = factor(param, levels = param_order)) %>%
    arrange(param, trans_idx)

  gt <- as.matrix(ground_truth)
  j  <- match(df$param, colnames(gt))
  ok <- !is.na(j) & df$trans_idx >= 1 & df$trans_idx <= nrow(gt)

  gt_vals <- rep(NA_real_, nrow(df))
  gt_vals[ok] <- gt[cbind(df$trans_idx[ok], j[ok])]

  df$gt <- gt_vals
  df$rel_bias <- if_else(gt_vals != 0, (df$est- gt_vals) / abs(gt_vals), NA_real_)
  df$coverage <- as.integer( gt_vals >= df$lwr & gt_vals <= df$upr)

  df$typeI    <- ifelse( gt_vals == 0,
                         as.integer(df$lwr > 0 | df$upr < 0),
                         NA_integer_)

  df$typeII <- ifelse(abs(gt_vals) > 0,
                      as.integer(df$lwr < 0  &  df$upr > 0), NA_integer_)

  df


}



# Bias-centered Rubin CIs
Ic_computation <- function(avg_parameters, all_fits,
                           bias_hat   = NULL,
                           alpha      = 0.05,
                           use_t      = TRUE) {
  stopifnot(length(all_fits) >= 2)
  m <- length(all_fits)
  p <- ncol(avg_parameters)
  par_names <- colnames(avg_parameters)

  # -------- 1) Collect aligned coef & vcov for each imputation ----------
  param_matrix <- vector("list", m)
  U_list       <- vector("list", m)

  for (i in seq_len(m)) {
    coefs_i <- lapply(all_fits[[i]], stats::coef)
    par_i <- sapply(1:3, function(k) {
      v <- coefs_i[[k]]
      z <- setNames(rep(NA_real_, p), par_names)
      if (!is.null(v)) {
        nm <- intersect(names(v), par_names)
        z[nm] <- v[nm]
      }
      z
    })
    colnames(par_i) <- paste0("model", 1:3)
    param_matrix[[i]] <- par_i

    vcovs_i <- lapply(all_fits[[i]], stats::vcov)
    U_list[[i]] <- lapply(vcovs_i, function(V) {
      V2 <- matrix(0, p, p, dimnames = list(par_names, par_names))
      if (!is.null(V)) {
        rn <- intersect(rownames(V), par_names)
        cn <- intersect(colnames(V), par_names)
        V2[rn, cn] <- V[rn, cn, drop = FALSE]
      }
      V2
    })
  }

  # -------- 2) Within-imputation variance U_bar ----------
  U_bar <- lapply(1:3, function(k) Reduce("+", lapply(U_list, `[[`, k)) / m)

  # -------- 3) Between-imputation variance B ----------
  B <- vector("list", 3)
  for (k in 1:3) {
    Theta_k <- sapply(param_matrix, function(M) M[, k])     # p x m
    center  <- avg_parameters[k, ]                          # pooled mean across imputations
    Ck      <- sweep(Theta_k, 1, center, "-")               # p x m
    B[[k]]  <- (Ck %*% t(Ck)) / (m - 1)                     # p x p
  }

  # -------- 4) Rubin total variance T = U_bar + (1 + 1/m) * B ----------
  Tm <- lapply(1:3, function(k) U_bar[[k]] + (1 + 1/m) * B[[k]])

  # -------- 5) Quantiles & df (optional t; otherwise z) ----------
  get_crit <- function(k) {
    if (!use_t) {
      return(rep(stats::qnorm(1 - alpha/2), p))
    }
    Ub <- diag(U_bar[[k]]); Ub[Ub <= 0] <- NA_real_
    Bb <- pmax(diag(B[[k]]), 0)
    r  <- (1 + 1/m) * Bb / Ub; r[!is.finite(r)] <- 0
    nu_rubin <- (m - 1) * (1 + 1/r)^2
    stats::qt(1 - alpha/2, df = nu_rubin)
  }

  crit <- lapply(1:3, get_crit)
  SE   <- lapply(1:3, function(k) sqrt(diag(Tm[[k]])))

  # -------- 6) Center: bias-corrected if 'bias_hat' provided ----------
  if (is.null(bias_hat)) {
    theta_center <- avg_parameters
  } else {

    stopifnot(nrow(bias_hat) == 3, ncol(bias_hat) == p)
    colnames(bias_hat) <- par_names
    theta_center <- avg_parameters - bias_hat
  }

  # -------- 7) Build CIs ----------
  CI <- vector("list", 3); names(CI) <- paste0("model", 1:3)
  for (k in 1:3) {
    lower <- theta_center[k, ] - crit[[k]] * SE[[k]]
    upper <- theta_center[k, ] + crit[[k]] * SE[[k]]
    CI[[k]] <- cbind(`L95%` = lower, `U95%` = upper)
    rownames(CI[[k]]) <- par_names
  }
  return(CI)
}


ci_dispatch$mipd_iter <- function(obj, bias_hat =NULL) {


  est <- averaging_params(obj)
  ci  <- Ic_computation(est, obj)

  df <- map_dfr(seq_along(ci), function(k) {
    ci_k   <- as.matrix(ci[[k]])
    pars   <- rownames(ci_k)
    est_k  <- est[k,]
    est_k  <- setNames(as.numeric(est_k), names(est_k))

    tibble(
      trans_idx = k,
      param     = pars,
      est       = unname(est_k[pars]),
      lwr       = ci_k[, 1],
      upr       = ci_k[, 2]
    )
  })

  param_order <- c("shape", "rate", "cov")
  df <- df %>%
    mutate(param = factor(param, levels = param_order)) %>%
    arrange(param, trans_idx)

  gt <- as.matrix(ground_truth)
  j  <- match(df$param, colnames(gt))
  ok <- !is.na(j) & df$trans_idx >= 1 & df$trans_idx <= nrow(gt)

  gt_vals <- rep(NA_real_, nrow(df))
  gt_vals[ok] <- gt[cbind(df$trans_idx[ok], j[ok])]

  df$gt <- gt_vals
  df$rel_bias <- if_else(gt_vals != 0, (df$est- gt_vals) / abs(gt_vals), NA_real_)
  df$coverage <- as.integer(!is.na(gt_vals) & gt_vals >= df$lwr & gt_vals <= df$upr)

  df$typeI    <- ifelse( gt_vals == 0,
                         as.integer(df$lwr > 0 | df$upr < 0),
                         NA_integer_)

  df$typeII <- ifelse(abs(gt_vals) > 0,
                      as.integer(df$lwr < 0  &  df$upr > 0), NA_integer_)

  df


}



averaging_params <- function(all_fits) {
  param_names <- names(stats::coef(all_fits[[1]][[1]]))
  out <- matrix(0, nrow = length(all_fits[[1]]), ncol = length(param_names))
  colnames(out) <- param_names
  for (i in seq_along(all_fits[[1]])) {
    for (p in seq_along(param_names)) {
      out[i, p] <- mean(sapply(all_fits, function(fit) stats::coef(fit[[i]])[p]))
    }
  }
  out
}

Ic_computation_cate <- function(all_fits, alpha = 0.05, use_t = TRUE) {
  stopifnot(length(all_fits) >= 2)
  avg_parameters <- averaging_params(all_fits)
  m <- length(all_fits)
  #p <- length(avg_parameters)
  #print(p)
  p <- ncol(avg_parameters)
  par_names <- colnames(avg_parameters)
  param_matrix <- vector("list", m)
  U_list       <- vector("list", m)

  for (i in seq_len(m)) {
    coefs_i <- lapply(all_fits[[i]], stats::coef)
    par_i <- sapply(1:3, function(k) {
      v <- coefs_i[[k]]
      z <- setNames(rep(NA_real_, p), par_names)
      if (!is.null(v)) {
        nm <- intersect(names(v), par_names)
        z[nm] <- v[nm]
      }
      z
    })
    colnames(par_i) <- paste0("model", 1:3)
    param_matrix[[i]] <- par_i

    vcovs_i <- lapply(all_fits[[i]], stats::vcov)
    U_list[[i]] <- lapply(vcovs_i, function(V) {
      V2 <- matrix(0, p, p, dimnames = list(par_names, par_names))
      if (!is.null(V)) {
        rn <- intersect(rownames(V), par_names)
        cn <- intersect(colnames(V), par_names)
        V2[rn, cn] <- V[rn, cn, drop = FALSE]
      }
      V2
    })
  }

  U_bar <- lapply(1:3, function(k) Reduce("+", lapply(U_list, `[[`, k)) / m)

  B <- vector("list", 3)
  for (k in 1:3) {
    Theta_k <- sapply(param_matrix, function(M) M[, k])   # p x m
    center  <- avg_parameters[k, ]
    Ck      <- sweep(Theta_k, 1, center, "-")             # p x m
    B[[k]]  <- (Ck %*% t(Ck)) / (m - 1)                   # p x p
  }

  Tm <- lapply(1:3, function(k) U_bar[[k]] + (1 + 1/m) * B[[k]])

  get_crit <- function(k) {
    if (!use_t) return(rep(stats::qnorm(1 - alpha/2), p))
    Ub <- diag(U_bar[[k]]); Ub[Ub <= 0] <- NA_real_
    Bb <- pmax(diag(B[[k]]), 0)
    r  <- (1 + 1/m) * Bb / Ub; r[!is.finite(r)] <- 0
    nu_rubin <- (m - 1) * (1 + 1/r)^2
    stats::qt(1 - alpha/2, df = nu_rubin)
  }
  crit <- lapply(1:3, get_crit)
  SE   <- lapply(1:3, function(k) sqrt(diag(Tm[[k]])))

  CI <- vector("list", 3); names(CI) <- paste0("model", 1:3)
  for (k in 1:3) {
    lower <- avg_parameters[k, ] - crit[[k]] * SE[[k]]
    upper <- avg_parameters[k, ] + crit[[k]] * SE[[k]]
    CI[[k]] <- cbind(est=avg_parameters[k,],`L95%` = lower, `U95%` = upper)
    rownames(CI[[k]]) <- colnames(avg_parameters)
  }
  CI
}

Ic_computation <- function(avg_parameters, all_fits, alpha = 0.05, use_t = TRUE) {
  stopifnot(length(all_fits) >= 2)
  m <- length(all_fits)
  p <- length(avg_parameters)
  #p <- ncol(avg_parameters)
  par_names <- colnames(avg_parameters)
  param_matrix <- vector("list", m)
  U_list       <- vector("list", m)

  for (i in seq_len(m)) {
    coefs_i <- lapply(all_fits[[i]], stats::coef)
    par_i <- sapply(1:3, function(k) {
      v <- coefs_i[[k]]
      z <- setNames(rep(NA_real_, p), par_names)
      if (!is.null(v)) {
        nm <- intersect(names(v), par_names)
        z[nm] <- v[nm]
      }
      z
    })
    colnames(par_i) <- paste0("model", 1:3)
    param_matrix[[i]] <- par_i

    vcovs_i <- lapply(all_fits[[i]], stats::vcov)
    U_list[[i]] <- lapply(vcovs_i, function(V) {
      V2 <- matrix(0, p, p, dimnames = list(par_names, par_names))
      if (!is.null(V)) {
        rn <- intersect(rownames(V), par_names)
        cn <- intersect(colnames(V), par_names)
        V2[rn, cn] <- V[rn, cn, drop = FALSE]
      }
      V2
    })
  }

  U_bar <- lapply(1:3, function(k) Reduce("+", lapply(U_list, `[[`, k)) / m)

  B <- vector("list", 3)
  for (k in 1:3) {
    Theta_k <- sapply(param_matrix, function(M) M[, k])   # p x m
    center  <- avg_parameters[k, ]
    Ck      <- sweep(Theta_k, 1, center, "-")             # p x m
    B[[k]]  <- (Ck %*% t(Ck)) / (m - 1)                   # p x p
  }

  Tm <- lapply(1:3, function(k) U_bar[[k]] + (1 + 1/m) * B[[k]])

  get_crit <- function(k) {
    if (!use_t) return(rep(stats::qnorm(1 - alpha/2), p))
    Ub <- diag(U_bar[[k]]); Ub[Ub <= 0] <- NA_real_
    Bb <- pmax(diag(B[[k]]), 0)
    r  <- (1 + 1/m) * Bb / Ub; r[!is.finite(r)] <- 0
    nu_rubin <- (m - 1) * (1 + 1/r)^2
    stats::qt(1 - alpha/2, df = nu_rubin)
  }
  crit <- lapply(1:3, get_crit)
  SE   <- lapply(1:3, function(k) sqrt(diag(Tm[[k]])))

  CI <- vector("list", 3); names(CI) <- paste0("model", 1:3)
  for (k in 1:3) {
    lower <- avg_parameters[k, ] - crit[[k]] * SE[[k]]
    upper <- avg_parameters[k, ] + crit[[k]] * SE[[k]]
    CI[[k]] <- cbind(`L95%` = lower, `U95%` = upper)
    rownames(CI[[k]]) <- colnames(avg_parameters)
  }
  CI
}

ci_teach_stud <- function(all_fits, ground_truth, alpha = 0.05, use_t = TRUE, param_order = NULL) {
  est <- averaging_params(all_fits)
  ci  <- Ic_computation(est, all_fits, alpha = alpha, use_t = use_t)
  df <- purrr::map_dfr(seq_along(ci), function(k) {
    ci_k   <- as.matrix(ci[[k]])
    pars   <- rownames(ci_k)
    est_k  <- setNames(as.numeric(est[k, ]), colnames(est))
    tibble::tibble(
      trans_idx = k,
      param     = pars,
      est       = unname(est_k[pars]),
      lwr       = ci_k[, 1],
      upr       = ci_k[, 2]
    )
  })
  if (!is.null(param_order)) {
    df <- df %>%
      dplyr::mutate(param = factor(param, levels = param_order)) %>%
      dplyr::arrange(param, trans_idx)
  }
  gt <- as.matrix(ground_truth)
  j  <- match(as.character(df$param), colnames(gt))
  ok <- !is.na(j) & df$trans_idx >= 1 & df$trans_idx <= nrow(gt)
  gt_vals <- rep(NA_real_, nrow(df))
  gt_vals[ok] <- gt[cbind(df$trans_idx[ok], j[ok])]
  df$gt <- gt_vals
  df$rel_bias <- dplyr::if_else(!is.na(gt_vals) & gt_vals != 0, (df$est - gt_vals) / abs(gt_vals), NA_real_)
  df$coverage <- as.integer(!is.na(gt_vals) & gt_vals >= df$lwr & gt_vals <= df$upr)
  df$typeI    <- dplyr::if_else(!is.na(gt_vals) & gt_vals == 0,
                                as.integer(df$lwr > 0 | df$upr < 0),
                                NA_integer_)
  df$typeII   <- dplyr::if_else(!is.na(gt_vals) & gt_vals != 0,
                                as.integer(df$lwr < 0 & df$upr > 0),
                                NA_integer_)
  df
}


compute_ci_summaries <- function(
    sim_root,
    out_root,
    scheme_ids,
    term_rate,
    model_tags,
    ci_dispatch,
    gt_dir            = NULL,
    seeds             = NULL,
    filename          = "fits.RData",
    alpha             = 0.05,
    use_t             = TRUE,
    auto_install_pkgs = FALSE,
    attach_pkgs       = TRUE,
    overwrite         = TRUE,
    verbose           = TRUE
) {
  # ----------------- dependency handling -----------------
  .req <- c("tibble","dplyr","purrr","stringr")
  .ensure <- function(pkgs, auto_install = FALSE, attach = TRUE) {
    miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
    if (length(miss)) {
      if (!auto_install) {
        stop("Missing packages: ", paste(miss, collapse = ", "),
             ". Install them or set auto_install_pkgs=TRUE.")
      }
      message("Installing missing packages: ", paste(miss, collapse = ", "))
      install.packages(miss, repos = "https://cloud.r-project.org")
    }
    if (attach) {
      suppressPackageStartupMessages({
        lapply(pkgs, function(p) do.call("library", list(p, character.only = TRUE)))
      })
    }
  }
  .ensure(.req, auto_install_pkgs, attach_pkgs)

  `%>%` <- dplyr::`%>%`

  # ----------------- tiny helpers -----------------
  parse_seed <- function(path_or_name) {
    b <- basename(path_or_name)
    m <- stringr::str_match(b, "(?i)seed[_-]?([0-9]+)")
    if (!is.na(m[1,2])) return(as.integer(m[1,2]))
    NA_integer_
  }

  .parse_seed_dir <- function(x) {
    m <- regexec("^seed_([0-9]+)$", x)
    hit <- regmatches(x, m)[[1]]
    if (length(hit) == 2) as.integer(hit[2]) else NA_integer_
  }

  list_model_files <- function(sim_root, scheme_id, term_rate, tag,
                               filename = "fits.RData", seeds = NULL) {
    stopifnot(is.character(sim_root), dir.exists(sim_root))
    base_dir <- file.path(sim_root, sprintf("cov_scheme_%d", scheme_id), term_rate, tag)
    if (!dir.exists(base_dir)) {
      if (isTRUE(verbose)) warning("Base dir not found: ", base_dir)
      return(tibble::tibble(filepath=character(), seed=integer(), model_sub=character()))
    }
    # only immediate seed_* subdirs
    sdirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)
    sdirs <- sdirs[grepl("/seed_[0-9]+$", sdirs)]
    if (!is.null(seeds)) {
      keep_seed <- vapply(basename(sdirs), .parse_seed_dir, integer(1))
      sdirs <- sdirs[keep_seed %in% seeds]
    }
    files <- file.path(sdirs, "results", filename)
    exists <- file.exists(files)
    if (!any(exists)) {
      return(tibble::tibble(filepath=character(), seed=integer(), model_sub=character()))
    }
    files <- files[exists]
    seeds_found <- vapply(basename(dirname(dirname(files))), .parse_seed_dir, integer(1))
    tibble::tibble(
      filepath = files,
      seed     = seeds_found,
      model_sub = tag
    ) %>%
      dplyr::mutate(seed = dplyr::if_else(is.na(seed), dplyr::row_number(), seed)) %>%
      dplyr::arrange(seed)
  }

  safe_ci_call <- function(tag, obj) {
    fn <- ci_dispatch[[tag]]
    if (!is.function(fn))
      stop(sprintf("No CI function provided for tag '%s'", tag))
    out <- fn(obj)
    # must at least have these
    req <- c("param","est","lwr","upr")
    if (!all(req %in% names(out)))
      stop(sprintf("CI fn for '%s' must return columns: %s",
                   tag, paste(req, collapse=", ")))
    out <- tibble::as_tibble(out)
    if (!("trans_idx" %in% names(out))) out$trans_idx <- NA_integer_
    out$param <- as.character(out$param)
    out
  }



  summarize_across_seeds <- function(per_seed) {
    per_seed %>%
      dplyr::group_by(scheme_id, model, param, trans_idx) %>%
      dplyr::summarise(
        n_seeds       = dplyr::n_distinct(seed),
        est_mean      = mean(est, na.rm = TRUE),
        est_sd        = stats::sd(est, na.rm = TRUE),
        rel_bias      = mean(rel_bias, na.rm = TRUE),
        coverage_rate = mean(coverage, na.rm = TRUE),
        type_I_rate   = mean(typeI, na.rm = TRUE),
        type_II_rate  = mean(typeII, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::arrange(param, trans_idx, model, scheme_id)
  }

  # ----------------- output root -----------------
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  # ----------------- main loop -----------------
  all_outputs <- list()

  for (scheme_id in scheme_ids) {
    if (isTRUE(verbose)) cat("\n=== Scheme", scheme_id, "===\n")
    # ground truth for this scheme (optional)
    gt_path <- if (!is.null(gt_dir)) {
      file.path(gt_dir, term_rate, sprintf("params_scheme_%02d.rds", scheme_id))
    } else NULL
    ground_truth <- if (!is.null(gt_path) && file.exists(gt_path)) {
      readRDS(gt_path)
    } else {
      if (!is.null(gt_dir) && isTRUE(verbose))
        warning(sprintf("[scheme %d] ground truth not found at %s", scheme_id, gt_path))
      NULL
    }

    scheme_out_dir <- file.path(out_root, sprintf("cov_scheme_%d", scheme_id), term_rate)
    dir.create(scheme_out_dir, recursive = TRUE, showWarnings = FALSE)

    per_tag_results <- list()

    for (tag in model_tags) {
      if (isTRUE(verbose)) cat(" • Model tag:", tag, "\n")

      out_base <- file.path(out_root, sprintf("cov_scheme_%d", scheme_id), term_rate, tag)
      ic_path  <- file.path(out_base, "ic_summary.rds")
      if (file.exists(ic_path) && !overwrite) {
        if (isTRUE(verbose)) cat("   - summary exists, skipping (set overwrite=TRUE):", ic_path, "\n")
        per_tag_results[[tag]] <- readRDS(ic_path)
        next
      }

      files_df <- list_model_files(sim_root, scheme_id, term_rate, tag, filename = filename, seeds = seeds)
      if (nrow(files_df) == 0) {
        warning(sprintf("[scheme %d | %s] no files found", scheme_id, tag))
        next
      }

      # per-seed extraction
      per_seed <- purrr::pmap_dfr(
        files_df,
        function(filepath, seed, model_sub = NA_character_, kind = NULL) {
          # Load model
          obj_name <- load(filepath)
          obj <- get(obj_name)

          # Directly call CI computation
          ci <- try(ci_teach_stud(obj, ground_truth), silent = TRUE)

          if (inherits(ci, "try-error") || nrow(ci) == 0) {
            warning(sprintf("[scheme %d | %s] seed %s: CI computation failed", scheme_id, tag, seed))
            return(tibble::tibble())
          }


          # Tag info
          reported_model <- if (!is.na(model_sub)) model_sub else tag
          dplyr::mutate(ci, scheme_id = scheme_id, model = reported_model, seed = seed)
        }
      )


      if (nrow(per_seed) == 0) {
        warning(sprintf("[scheme %d | %s] nothing extracted", scheme_id, tag))
        next
      }

      ic_summary <- summarize_across_seeds(per_seed)
      dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
      saveRDS(ic_summary, file = ic_path)

      per_tag_results[[tag]] <- ic_summary
    } # /for tag

    all_outputs[[as.character(scheme_id)]] <- per_tag_results
  } # /for scheme


}

ci_dispatch <- list(
  dropout_adaptive_tau = ci_teach_stud
)
