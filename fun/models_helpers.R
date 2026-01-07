# cox

prepare_cox <- function(df){

  scheme_visits <- df %>% dplyr::arrange(patient_id, visits)
  n_pats <- dplyr::n_distinct(scheme_visits$patient_id)


  scheme_data <- scheme_visits %>%
    dplyr::group_by(patient_id) %>%
    dplyr::mutate(onset = as.integer(any(onset == 1))) %>%
    dplyr::ungroup()


  row_id <- scheme_visits %>%
    dplyr::count(patient_id, name = "nrows") %>%
    dplyr::arrange(patient_id)

  scheme_visits$patient_id <- rep(seq_len(n_pats), times = as.numeric(row_id$nrows))
  scheme_data$patient_id   <- rep(seq_len(n_pats), times = as.numeric(row_id$nrows))


  scheme_data <- scheme_data %>%
    dplyr::group_by(patient_id) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    dplyr::select(-visits)

  as.data.frame(scheme_data)
}

fit_cox <- function(df,covariate_names="cov") {
  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()),
                           names = c("Disease-free", "Disease", "Death"))

  # covariate_names <- setdiff(names(df),
  #                            c("patient_id","age","onset","onset_age",
  #                              "dead","death_time","last_bfo","visits"))

  data_long <- mstate::msprep(
    data   = df,
    trans  = tmat,
    time   = c(NA, "onset_age", "death_time"),
    status = c(NA, "onset", "dead"),
    keep   = c("age", covariate_names),
    id     = "patient_id"
  )

  data_long$Tstart[data_long$trans < 3] <- data_long$Tstart[data_long$trans < 3] + data_long$age[data_long$trans < 3]
  data_long$time <- data_long$Tstop - data_long$Tstart

  data_long <- mstate::expand.covs(data_long, covariate_names)

  expanded_covariates <- setdiff(names(data_long),
                                 c("patient_id","from","to","trans","Tstart","Tstop","time","status","age","cov"))

  fml <- stats::as.formula(
    paste("survival::Surv(Tstart, Tstop, status) ~",
          paste(expanded_covariates, collapse = " + "),
          "+ strata(trans)")
  )

  # Return model (not just coefs)
  survival::coxph(fml, data = data_long, method = "breslow")
}

fit_one_seed_cox <- function(seed,path,out_dir,model_tag) {
  sim_path <- file.path(path,
                        sprintf("simulation_ready_%03d.rds", seed))
  print(sim_path)
  if (!file.exists(sim_path)) {
    stop(sprintf("Missing file for seed %d -> skip", seed))
    return(NULL)
  }

  d <- readRDS(sim_path)
  panel_data <- d[[2]]
  temp <- prepare_cox(panel_data)

  # fit <-
  #   tryCatch(
  #   fit_cox(temp)
  #   , error = function(e) NULL)

  fit <- fit_cox(temp)

  coefs <- stats::coef(fit)


  out_dir <- file.path(out_dir,
                       model_tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


  saveRDS(fit,   file.path(out_dir, sprintf("seed_%03d_model_cox.rds", seed)))


  coefs
}

# flexsurv

prep_long <- function(df, tmat) {
  dl <- mstate::msprep(
    data   = df,
    trans  = tmat,
    time   = c(NA, "onset_age", "death_time"),
    status = c(NA, "onset", "dead"),
    keep   = c("age","cov"),
    id     = "patient_id"
  )
  dl$Tstart[dl$trans < 3] <- dl$Tstart[dl$trans < 3] + dl$age[dl$trans < 3]
  dl %>% mutate(time = Tstop - Tstart)
}

fit_one_seed_gomp <- function(seed,path,out_dir,model_tag) {
  sim_path <- file.path(path,
                         sprintf("simulation_ready_%03d.rds", seed))

  if (!file.exists(sim_path)) {
    message(sprintf("Missing file for seed %d -> skip", seed))
    return(NULL)
  }

  data <- readRDS(sim_path)
  real_data  <- data[[1]]
  panel_data <- data[[2]]

  # collapse panel to one row per patient
  panel_collapsed <- panel_data %>%
    group_by(patient_id) %>%
    mutate(onset = as.integer(any(onset == 1))) %>%
    slice(1) %>%
    ungroup() %>%
    select(-visits)

  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()),
                           names = c("Dementia-free","Dementia","Death"))

  data_long_real  <- prep_long(real_data,  tmat)
  data_long_panel <- prep_long(panel_collapsed, tmat)

  n_trans <- max(tmat, na.rm = TRUE)


  fits_real  <- vector("list", n_trans)
  fits_panel <- vector("list", n_trans)
  for (i in 1:n_trans) {
    sub_r <- subset(data_long_real,  trans == i)
    sub_p <- subset(data_long_panel, trans == i)

    fits_real[[i]] <- if (nrow(sub_r) == 0 || sum(sub_r$status) == 0) NA else
      tryCatch(flexsurvreg(Surv(Tstart, Tstop, status) ~ cov, data = sub_r, dist = "gompertz"),
               error = function(e) NA)

    fits_panel[[i]] <- if (nrow(sub_p) == 0 || sum(sub_p$status) == 0) NA else
      tryCatch(flexsurvreg(Surv(Tstart, Tstop, status) ~ cov, data = sub_p, dist = "gompertz"),
               error = function(e) NA)
  }


  get_param_names <- function(fits) {
    for (k in seq_along(fits)) {
      if (is.list(fits[[k]]) && !is.na(fits[[k]])[1]) return(names(fits[[k]]$coefficients))
    }

    return(c("shape","rate","cov"))
  }
  param_names <- get_param_names(fits_panel)


  to_mat <- function(fits, param_names) {
    M <- matrix(NA_real_, nrow = length(fits), ncol = length(param_names))
    colnames(M) <- param_names
    for (i in seq_along(fits)) {
      fi <- fits[[i]]
      if (is.list(fi) && !is.na(fi)[1]) {
        co <- fi$coefficients
        nm <- intersect(names(co), param_names)
        M[i, nm] <- unname(co[nm])
      }
    }
    rownames(M) <- seq_len(nrow(M))
    M
  }

  out_base <- file.path(out_dir,
                        model_tag)
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  saveRDS(fits_real,  file.path(out_base, sprintf("seed_%03d_model_flex_real.rds",  seed)))
  saveRDS(fits_panel,  file.path(out_base, sprintf("seed_%03d_model_flex_panel.rds",  seed)))

  params_real  <- to_mat(fits_real,  param_names)
  params_panel <- to_mat(fits_panel, param_names)

  list(panel = params_panel,real=params_real)
}

# msm / msm_age

prepare_msm <- function(df){
  df$state <- 1
  df <- df %>%
    arrange(patient_id, visits) %>%
    group_by(patient_id) %>%
    do({
      if (nrow(.) == 1) {
        duplicated_row <- bind_rows(., .)
        duplicated_row$age[2]    <- duplicated_row$death_time[1]
        duplicated_row$visits[2] <- 2
        duplicated_row
      } else {
        .
      }
    }) %>%
    ungroup()

  df <- df %>%
    group_by(patient_id) %>%
    mutate(
      state = ifelse(onset == 1, 2, state),
      state = ifelse(row_number() == n() & dead == 1, 3, state),
      state = ifelse(row_number() == n() & dead == 0 & onset == 0, 99, state)
    ) %>%
    ungroup()

  df
}

prepare_msm2 <- function(df){
  panel_data_death <- df%>% filter(dead==1 ) %>% group_by(patient_id) %>% mutate(age=death_time) %>% distinct(patient_id,.keep_all = T)
  df2 <- rbind(df,panel_data_death) %>% group_by(patient_id) %>%  mutate(state=case_when(max(onset)==1 & onset_age<=age& age!=death_time~2 ,
                                                              death_time==age & dead==1 ~3,
        TRUE~1)) %>% arrange(patient_id)
  return(df2)
}

fit_msm <- function(df) {
  Q <- rbind(c(0, 1, 1),
             c(0, 0, 1),
             c(0, 0, 0))

  fit <- tryCatch(
    msm(state ~ age,
        subject    = patient_id,
        data       = df,
        qmatrix    = Q,
        covariates = ~ cov,
        censor     = 99,
        control    = list(fnscale = 1000, maxit = 1000),
        deathexact = TRUE),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 2,
                  dimnames = list(c("1","2","3"),
                                  c("rate","cov")))
    return(list(model = NULL, params = mat))
  }

  v <- fit$estimates
  if (length(v) < 6) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 2,  byrow = FALSE,
                  dimnames = list(c("1","2","3"),
                                  c("rate","cov")))
  } else {
    mat <- matrix(v[1:6], nrow = 3, byrow = FALSE,
                  dimnames = list(c("1","2","3"),
                                  c("rate","cov")))
  }

  list(model = fit, params = mat)
}

fit_msm_age <- function(df,seed,out_dir) {

  msm_root <-  file.path(gsub("_age","",out_dir),
                          sprintf("seed_%03d_model_msm.rds", seed))
  model.msm <- readRDS(msm_root)
  #initial_guess_age <- qmatrix.msm(model.msm)$estimates
  Q <- rbind(c(0, 1, 1),
             c(0, 0, 1),
             c(0, 0, 0))


  fit <- tryCatch(
    msm(state ~ age,
        subject = patient_id,
        data = df,
        qmatrix = Q,
        gen.inits = F,
        covariates = ~ cov + age ,
        control    = list(fnscale = 1000, maxit = 1000),
        censor = 99,
        deathexact = TRUE),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 3,
                  dimnames = list(c("1","2","3"),
                                  c("rate","cov", "age")))
    return(list(model = NULL, params = mat))
  }

  v <- fit$estimates
  if (length(v) < 9) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 3, byrow = F,
                  dimnames = list(c("1","2","3"),
                                  c("rate","cov", "age")))
  } else {
    mat <- matrix(v[1:9], nrow = 3, byrow = F,
                  dimnames = list(c("1","2","3"),
                                  c("rate","cov", "age")))
  }

  list(model = fit, params = mat)
}

fit_one_seed_msm <- function(seed, age= FALSE,path,out_dir) {
  sim_path <- file.path(path,
                        sprintf("simulation_ready_%03d.rds", seed)
  )
  if (!file.exists(sim_path)) {
    message(sprintf("Missing file for seed %d -> skip", seed))
    return(NULL)
  }

  d <- readRDS(sim_path)
  panel_data <- d[[2]]
  temp_panel <- prepare_msm2(panel_data)
  out_base <- out_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (age == FALSE){
    res <- fit_msm(temp_panel)
    saveRDS(res$model,  file.path(out_base, sprintf("seed_%03d_model_msm.rds",  seed)))
  }else{
    res <- fit_msm_age(temp_panel, seed,out_base)
    saveRDS(res$model,  file.path(out_base, sprintf("seed_%03d_model_msm_age.rds",  seed)))
  }

  list( params = res$params)
}

# nhm

fit_nhm <- function(df) {
  tmat_1 <- rbind(c(0,1,2),c(0,0,3),rep(0,3))
  tmat_2 <- rbind(c(0,4,5),c(0,0,6),rep(0,3))
  tmat_3 <- rbind(c(0,7,8),c(0,0,9),rep(0,3))

  df$patient_id <- as.factor(df$patient_id)
  df %<>% group_by(patient_id) %>% filter(n()>1 )
  df$age <- round(df$age,3)
  df=as.data.frame(df)

  find_splits <- function(age) {
    quantiles <- quantile(age, probs = seq(0, 1, 0.02))
    #quantiles <-seq(min(age),max(age),by=0.001)
    return(quantiles[-c(1, length(quantiles))])
  }

  split_points <- find_splits(df$age)


  fit <- tryCatch({
    object_nhm <- model.nhm(
      state ~ age,
      subject       = patient_id,
      data          = df,
      trans         = tmat_1,
      nonh          = tmat_1,
      type          = "gompertz",
      covariates    = c("cov"),
      covm          = list(cov = tmat_1),
      # censor.states = c(1, 2),
      #  censor        = 99,
      death         = TRUE,
      death.states  = c(3)
    )

    model_nhm <- nhm(
      object_nhm,
      gen_inits   = TRUE,
      score_test  = FALSE,
      control     = nhm.control(ncores = 10,splits = split_points, obsinfo = FALSE)
    )

    model_nhm
  },
  error = function(e) NULL
  )


  if (is.null(fit)) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 3,
                  dimnames = list(c("1","2","3"),
                                  c("rate", "shape", "cov")))
    return(list(model = NULL, params = mat))
  }

  v <- fit$par
  if (length(v) < 9) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 3,
                  dimnames = list(c("1","2","3"),
                                  c("rate", "shape", "cov")))
  } else {
    mat <- matrix(v[1:9], nrow = 3, ncol = 3,
                  dimnames = list(c("1","2","3"),
                                  c("rate", "shape", "cov")))
  }

  list(model = fit, params = mat)
}

fit_one_seed_nhm <- function(seed, age= FALSE,path,out_dir) {
  sim_path <- file.path(path,
                        sprintf("simulation_ready_%03d.rds", seed))
  if (!file.exists(sim_path)) {
    message(sprintf("Missing file for seed %d -> skip", seed))
    return(NULL)
  }

  d <- readRDS(sim_path)
  panel_data <- d[[2]]
  temp_panel <- prepare_msm2(panel_data)

  res <- fit_nhm(temp_panel)
  out_base <- out_dir
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  saveRDS(res$model,  file.path(out_base, sprintf("seed_%03d_model_nhm.rds",  seed)))


  list( params = res$params, model = res$model)
}

# mipd

bias_bootstrap <- function(
    scheme_data, disease_status, disease_age,
    cov_vector, clock_assumption, distribution,
    averaging_params,
    B = 100, m_boot = 5L, seed = 1L
) {
  set.seed(seed)

  n_patients <- nrow(scheme_data)
  m <- ncol(disease_status)
  stopifnot(n_patients == nrow(disease_age), m == ncol(disease_age))

  use_j <- sort(sample.int(m, size = min(m_boot, m), replace = FALSE))

  boots <- vector("list", B)

  for (b in seq_len(B)) {

    samp_idx <- sample.int(n_patients, size = n_patients, replace = TRUE)

    fits_b <- lapply(use_j, function(j) {
      temp_b <- scheme_data[samp_idx, , drop = FALSE]

      temp_b$onset     <- disease_status[samp_idx, j]
      temp_b$onset_age <- disease_age[samp_idx, j]


      if ("patient_id" %in% names(temp_b)) {
        temp_b$patient_id <- seq_len(nrow(temp_b))
      }

      fit_model(temp_b, cov_vector, clock_assumption, distribution)
    })

    # Pool across the m_boot imputations for this bootstrap draw → 3 x p pooled estimate
    boots[[b]] <- averaging_params(fits_b)
  }

  boots
}

fit_model <- function(data, cov_vector, clock_assumption, distribution) {
  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()), names = c("Disease-free", "Disease", "Death"))


  covariate_names <- cov_vector


  data_long <- mstate::msprep(data = data, trans = tmat,
                              time = c(NA, "onset_age", "death_time"),
                              status = c(NA, "onset", "dead"),
                              keep = c("age", covariate_names),
                              id = "patient_id")

  data_long$Tstart[data_long$trans < 3] <- data_long$Tstart[data_long$trans < 3] + data_long$age[data_long$trans < 3]
  data_long$time <- data_long$Tstop - data_long$Tstart

  n_trans <- max(tmat, na.rm = TRUE)
  fits <- vector(mode = "list", length = n_trans)

  forward_formula <- stats::as.formula(paste("survival::Surv(Tstart, Tstop, status) ~", paste(covariate_names, collapse = " + ")))
  reset_formula <- stats::as.formula(paste("survival::Surv(time, status) ~", paste(covariate_names, collapse = " + ")))


  if (clock_assumption == "forward") {
    for (i in 1:3) {
      fits[[i]] <- flexsurv::flexsurvreg(forward_formula,
                                         data = subset(data_long, trans == i),
                                         dist = distribution)
    }
  } else if (clock_assumption == "mix") {
    for (i in 1:2) {
      fits[[i]] <- flexsurv::flexsurvreg(forward_formula,
                                         data = subset(data_long, trans == i),
                                         dist = distribution)
    }
    fits[[3]] <- flexsurv::flexsurvreg(reset_formula,
                                       data = subset(data_long, trans == 3),
                                       dist = distribution)
  }

  return(fits)
}

averaging_params <- function(all_fits) {

  # Extract parameter names from the first fitted model
  param_names <- names(stats::coef(all_fits[[1]][[1]]))
  n_params <- length(param_names)

  # Initialize matrix for averaged parameters
  averaged_params <- matrix(0, nrow = length(all_fits[[1]]), ncol = n_params)
  colnames(averaged_params) <- param_names

  # Loop over transitions and parameters
  for (i in seq_along(all_fits[[1]])) {
    for (p in seq_along(param_names)) {
      averaged_params[i, p] <- mean(sapply(all_fits, function(fit) stats::coef(fit[[i]])[p]))
    }
  }

  return(averaged_params)
}

run_imputation2 <- function(data, m, cov_vector,clock_assumption, distribution, inner_cores){

  # Step 1: Create base visit and per-patient data
  n_pats <- length(unique(data$patient_id))
  scheme_visits <- data %>% arrange(patient_id, visits)  # keep order stable
  scheme_data   <- scheme_visits
  original_index <- unique(scheme_data$patient_id)

  # ensure numeric vector
  scheme_data$last_bfo <- rep(NA_real_, nrow(scheme_data))

  for (i in original_index) {
    onset_idx <- which(scheme_data$onset[scheme_data$patient_id == i] == 1)[1]
    if (!is.na(onset_idx) && onset_idx > 1) {
      scheme_data$last_bfo[scheme_data$patient_id == i] <-
        scheme_data$age[scheme_data$patient_id == i][onset_idx - 1]
    }
  }

  # Mark patient as having onset = 1 if any visit has onset
  for (i in original_index) {
    if (any(scheme_data$onset[scheme_data$patient_id == i] == 1)) {
      scheme_data$onset[scheme_data$patient_id == i] <- 1
    }
  }

  scheme_data <- scheme_data[order(scheme_data$patient_id), ]

  # keep row counts aligned with patient ordering
  row_id <- scheme_visits %>%
    count(patient_id, name = "nrows") %>%
    arrange(patient_id)

  scheme_visits$patient_id <- rep(seq_len(n_pats), times = as.numeric(row_id$nrows))
  scheme_data$patient_id   <- rep(seq_len(n_pats), times = as.numeric(row_id$nrows))

  scheme_data <- scheme_data %>%
    group_by(patient_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(-visits)

  # Step 2: Define last_bfo correctly
  last_visits <- as.data.table(scheme_visits)[, .SD[which.max(visits)], by = patient_id]
  last_visits <- last_visits[, .(patient_id, visits, last_visit_age = age)]
  scheme_data <- merge(scheme_data, last_visits, by = "patient_id", all.x = TRUE, sort = FALSE)
  scheme_data$last_bfo <- ifelse(scheme_data$onset == 1, scheme_data$last_bfo, scheme_data$last_visit_age)
  scheme_data <- scheme_data %>% select(-last_visit_age)

  # Step 3: Define covariates and prepare output
  scheme_data <- as.data.frame(scheme_data)
  covariate_names <- cov_vector


  temp <- scheme_data

  # Define transition matrix and the number of transitions for multi-state model
  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()), names = c("Disease-free", "Disease", "Death"))
  n_trans <- max(tmat, na.rm = TRUE)  # Get the total number of transitions


  # Prepare the data for the multi-state model using msprep function
  data_long <- mstate::msprep(data = temp, trans = tmat,
                              time = c(NA, "onset_age", "death_time"),
                              status = c(NA, "onset", "dead"),
                              keep = c("age",covariate_names),
                              id = "patient_id")

  # Adjust the start times for transitions
  data_long$Tstart[data_long$trans < 3] <- data_long$Tstart[data_long$trans < 3] + data_long$age[data_long$trans < 3]
  data_long$time <- data_long$Tstop - data_long$Tstart

  # Expand covariates for the model
  data_long <- mstate::expand.covs(data_long, covariate_names)

  # Define covariates for Cox Proportional Hazards model
  expanded_covariates <- setdiff(names(data_long), c("patient_id", "from", "to", "trans", "Tstart", "Tstop", "time", "status", "age", covariate_names))

  # Create formula for Cox Proportional Hazards model
  formula_str <- paste("survival::Surv(Tstart, Tstop ,status) ~",
                       paste(expanded_covariates, collapse = " + "),
                       "+ strata(trans)")
  model_formula <- stats::as.formula(formula_str)

  # Fit Cox Proportional Hazards model using Breslow method
  model_cox <- survival::coxph(model_formula, data = data_long, method = "breslow")


  # Compute baseline hazards for each transition
  all_haz <- survival::basehaz(model_cox, center = FALSE)

  # Identify coefficient names and their indices
  coef_names <- names(model_cox$coefficients)


  # Identify coefficients for the different transitions
  to_onset_indices <- grep("\\.1$", coef_names)
  to_death_indices <- grep("\\.2$", coef_names)
  dis_to_death_indices <- grep("\\.3$", coef_names)

  # Extract corresponding coefficients for each transition
  onsetpar <- model_cox$coefficients[to_onset_indices]
  deathpar <- model_cox$coefficients[to_death_indices]
  deathpar_dis <- model_cox$coefficients[dis_to_death_indices]

  # Extract hazard values for each transition
  onsethaz <- all_haz[all_haz$strata == "trans=1", 1:2]
  deathhaz <- all_haz[all_haz$strata == "trans=2", 1:2]
  deathhaz_dis <- all_haz[all_haz$strata == "trans=3", 1:2]


  #Compute hazard
  compute_hazard <- function(df) {
    df <- rbind(data.frame(time = 0, hazard = 0), df)
    df$h <- c(diff(df$hazard) / diff(df$time), NA)
    return(df)
  }

  # Apply to all three hazard tables
  onsethaz <- compute_hazard(onsethaz)
  deathhaz <- compute_hazard(deathhaz)
  deathhaz_dis <- compute_hazard(deathhaz_dis)

  # Merge onset and death hazards
  combhaz <- merge(onsethaz, deathhaz, by = "time", all = TRUE, suffixes = c(".12", ".13"))

  # Ensure first values are not NA
  combhaz$hazard.12[1] <- ifelse(is.na(combhaz$hazard.12[1]), 0, combhaz$hazard.12[1])
  combhaz$hazard.13[1] <- ifelse(is.na(combhaz$hazard.13[1]), 0, combhaz$hazard.13[1])

  # Precompute time differences
  dt <- c(NA, diff(combhaz$time))

  # Fill NAs in hazard rates and recalculate cumulative hazards
  na_12 <- is.na(combhaz$h.12)
  na_13 <- is.na(combhaz$h.13)

  for (i in which(na_12)) {
    combhaz$h.12[i] <- combhaz$h.12[i - 1]
    combhaz$hazard.12[i] <- combhaz$hazard.12[i - 1] + combhaz$h.12[i - 1] * dt[i]
  }

  for (i in which(na_13)) {
    combhaz$h.13[i] <- combhaz$h.13[i - 1]
    combhaz$hazard.13[i] <- combhaz$hazard.13[i - 1] + combhaz$h.13[i - 1] * dt[i]
  }

  # Matrices to store imputed data
  disease_status <- matrix(rep(scheme_data$onset, m), nrow = nrow(scheme_data), ncol = m)
  disease_age <- matrix(0, nrow = nrow(scheme_data), ncol = m)
  set.seed(2)
  c1 <- matrix(stats::runif(nrow(scheme_data) * m), nrow = nrow(scheme_data), ncol = m)
  c2 <- matrix(stats::runif(nrow(scheme_data) * m), nrow = nrow(scheme_data), ncol = m)

  # Covariate matrix
  covar <- as.matrix(scheme_data[, covariate_names])

  # Linear predictor calculations
  lp.12 <- covar %*% as.vector(t(onsetpar))
  lp.13 <- covar %*% as.vector(t(deathpar))
  lp.23 <- covar %*% as.vector(t(deathpar_dis))

  n_patients <- nrow(scheme_data)
  results <- parallel::mclapply(1:n_patients, function(x)process_patient(x,combhaz$time, combhaz$hazard.12, combhaz$hazard.13, combhaz$h.12,
                                                                           lp.12, lp.13, c1, c2,
                                                                           scheme_visits[scheme_visits$patient_id == x, ],
                                                                           m, mean(scheme_data$onset)), mc.cores = inner_cores)


  disease_age <- do.call(rbind, lapply(results, `[[`, "age"))
  disease_status <- do.call(rbind, lapply(results, `[[`, "status"))

  all_fits <-  parallel::mclapply(1:m, function(j) {
    temp <- scheme_data
    disease_generated <- which(disease_status[, j] == 1)
    temp$onset[disease_generated] <- 1
    temp$onset_age[disease_generated] <- disease_age[disease_generated, j]
    fit_model(temp, cov_vector, clock_assumption, distribution)
  })



  # Compute averaged parameters
  averaged_params <- averaging_params(all_fits)

  return(list(all_fits, averaged_params))
}


fit_one_seed_mipd <- function(seed,cov_vector,path,out_dir) {
  sim_path <- file.path(path,
                        sprintf("simulation_ready_%03d.rds", seed))
  if (!file.exists(sim_path)) {
    message(sprintf("Missing file for seed %d -> skip", seed))
    return(NULL)
  }

  d <- readRDS(sim_path)
  panel_data <- d[[2]]

  m <- 30
  clock_assumption <- "forward"
  distribution <- "gompertz"
  inner_cores<-20

  res <- tryCatch(
    run_imputation2(panel_data, m, cov_vector, clock_assumption, distribution, inner_cores),
    error = function(e) NULL
  )


  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(res[[1]],  file.path(out_dir, sprintf("seed_%03d_model_mipd.rds",  seed)))


  list(params = res[[2]])
}


# mipd iter

run_imputation_iterative <- function(data, m, cov_vector, clock_assumption, distribution, inner_cores, max_iter, eps, boot = FALSE){

  # Step 1: Create base visit and per-patient data
  n_pats <- length(unique(data$patient_id))
  scheme_visits <- data %>% arrange(patient_id, visits)  # keep order stable
  scheme_data   <- scheme_visits
  original_index <- unique(scheme_data$patient_id)

  # ensure numeric vector
  scheme_data$last_bfo <- rep(NA_real_, nrow(scheme_data))

  for (i in original_index) {
    onset_idx <- which(scheme_data$onset[scheme_data$patient_id == i] == 1)[1]
    if (!is.na(onset_idx) && onset_idx > 1) {
      scheme_data$last_bfo[scheme_data$patient_id == i] <-
        scheme_data$age[scheme_data$patient_id == i][onset_idx - 1]
    }
  }

  # Mark patient as having onset = 1 if any visit has onset
  for (i in original_index) {
    if (any(scheme_data$onset[scheme_data$patient_id == i] == 1)) {
      scheme_data$onset[scheme_data$patient_id == i] <- 1
    }
  }

  scheme_data <- scheme_data[order(scheme_data$patient_id), ]

  # keep row counts aligned with patient ordering
  row_id <- scheme_visits %>%
    count(patient_id, name = "nrows") %>%
    arrange(patient_id)

  scheme_visits$patient_id <- rep(seq_len(n_pats), times = as.numeric(row_id$nrows))
  scheme_data$patient_id   <- rep(seq_len(n_pats), times = as.numeric(row_id$nrows))

  scheme_data <- scheme_data %>%
    group_by(patient_id) %>%
    slice(1) %>%
    ungroup() %>%
    select(-visits)

  # Step 2: Define last_bfo correctly
  last_visits <- as.data.table(scheme_visits)[, .SD[which.max(visits)], by = patient_id]
  last_visits <- last_visits[, .(patient_id, visits, last_visit_age = age)]
  scheme_data <- merge(scheme_data, last_visits, by = "patient_id", all.x = TRUE, sort = FALSE)
  scheme_data$last_bfo <- ifelse(scheme_data$onset == 1, scheme_data$last_bfo, scheme_data$last_visit_age)
  scheme_data <- scheme_data %>% select(-last_visit_age)

  # Step 3: Define covariates and prepare output
  scheme_data <- as.data.frame(scheme_data)
  covariate_names <- cov_vector


  temp <- scheme_data

  # Define transition matrix and the number of transitions for multi-state model
  tmat <- mstate::transMat(x = list(c(2, 3), c(3), c()), names = c("Disease-free", "Disease", "Death"))
  n_trans <- max(tmat, na.rm = TRUE)  # Get the total number of transitions


  # Prepare the data for the multi-state model using msprep function
  data_long <- mstate::msprep(data = temp, trans = tmat,
                              time = c(NA, "onset_age", "death_time"),
                              status = c(NA, "onset", "dead"),
                              keep = c("age",covariate_names),
                              id = "patient_id")

  # Adjust the start times for transitions
  data_long$Tstart[data_long$trans < 3] <- data_long$Tstart[data_long$trans < 3] + data_long$age[data_long$trans < 3]
  data_long$time <- data_long$Tstop - data_long$Tstart

  # Expand covariates for the model
  data_long <- mstate::expand.covs(data_long, covariate_names)

  # Define covariates for Cox Proportional Hazards model
  expanded_covariates <- setdiff(names(data_long), c("patient_id", "from", "to", "trans", "Tstart", "Tstop", "time", "status", "age", covariate_names))

  # Create formula for Cox Proportional Hazards model
  formula_str <- paste("survival::Surv(Tstart, Tstop ,status) ~",
                       paste(expanded_covariates, collapse = " + "),
                       "+ strata(trans)")
  model_formula <- stats::as.formula(formula_str)

  # Fit Cox Proportional Hazards model using Breslow method
  model_cox <- survival::coxph(model_formula, data = data_long, method = "breslow")


  # Compute baseline hazards for each transition
  all_haz <- survival::basehaz(model_cox, center = FALSE)

  # Identify coefficient names and their indices
  coef_names <- names(model_cox$coefficients)


  # Identify coefficients for the different transitions
  to_onset_indices <- grep("\\.1$", coef_names)
  to_death_indices <- grep("\\.2$", coef_names)
  dis_to_death_indices <- grep("\\.3$", coef_names)

  # Extract corresponding coefficients for each transition
  onsetpar <- model_cox$coefficients[to_onset_indices]
  deathpar <- model_cox$coefficients[to_death_indices]
  deathpar_dis <- model_cox$coefficients[dis_to_death_indices]

  # Extract hazard values for each transition
  onsethaz <- all_haz[all_haz$strata == "trans=1", 1:2]
  deathhaz <- all_haz[all_haz$strata == "trans=2", 1:2]
  deathhaz_dis <- all_haz[all_haz$strata == "trans=3", 1:2]


  #Compute hazard
  compute_hazard <- function(df) {
    df <- rbind(data.frame(time = 0, hazard = 0), df)
    df$h <- c(diff(df$hazard) / diff(df$time), NA)
    return(df)
  }

  # Apply to all three hazard tables
  onsethaz <- compute_hazard(onsethaz)
  deathhaz <- compute_hazard(deathhaz)
  deathhaz_dis <- compute_hazard(deathhaz_dis)

  # Merge onset and death hazards
  combhaz <- merge(onsethaz, deathhaz, by = "time", all = TRUE, suffixes = c(".12", ".13"))

  # Ensure first values are not NA
  combhaz$hazard.12[1] <- ifelse(is.na(combhaz$hazard.12[1]), 0, combhaz$hazard.12[1])
  combhaz$hazard.13[1] <- ifelse(is.na(combhaz$hazard.13[1]), 0, combhaz$hazard.13[1])

  # Precompute time differences
  dt <- c(NA, diff(combhaz$time))

  # Fill NAs in hazard rates and recalculate cumulative hazards
  na_12 <- is.na(combhaz$h.12)
  na_13 <- is.na(combhaz$h.13)

  for (i in which(na_12)) {
    combhaz$h.12[i] <- combhaz$h.12[i - 1]
    combhaz$hazard.12[i] <- combhaz$hazard.12[i - 1] + combhaz$h.12[i - 1] * dt[i]
  }

  for (i in which(na_13)) {
    combhaz$h.13[i] <- combhaz$h.13[i - 1]
    combhaz$hazard.13[i] <- combhaz$hazard.13[i - 1] + combhaz$h.13[i - 1] * dt[i]
  }

  # Matrices to store imputed data
  disease_status <- matrix(rep(scheme_data$onset, m), nrow = nrow(scheme_data), ncol = m)
  disease_age <- matrix(0, nrow = nrow(scheme_data), ncol = m)
  set.seed(2)
  c1 <- matrix(stats::runif(nrow(scheme_data) * m), nrow = nrow(scheme_data), ncol = m)
  c2 <- matrix(stats::runif(nrow(scheme_data) * m), nrow = nrow(scheme_data), ncol = m)

  # Covariate matrix
  covar <- as.matrix(scheme_data[, covariate_names])

  # Linear predictor calculations
  lp.12 <- covar %*% as.vector(t(onsetpar))
  lp.13 <- covar %*% as.vector(t(deathpar))
  lp.23 <- covar %*% as.vector(t(deathpar_dis))


  criteria <- Inf
  old <- list(onsetpar, deathpar, deathpar_dis)
  k <- 0
  while (criteria > eps && k < max_iter) {


    n_patients <- nrow(scheme_data)
    results <- parallel::mclapply(1:n_patients, function(x)process_patient(x,combhaz$time, combhaz$hazard.12, combhaz$hazard.13, combhaz$h.12,
                                                                           lp.12, lp.13, c1, c2,
                                                                           scheme_visits[scheme_visits$patient_id == x, ],
                                                                           m, mean(scheme_data$onset)), mc.cores = inner_cores)


    disease_age <- do.call(rbind, lapply(results, `[[`, "age"))
    disease_status <- do.call(rbind, lapply(results, `[[`, "status"))


    result_iteration <- parallel::mclapply(1:m, function(j) {
      temp <- scheme_data
      disease_generated <- which(disease_status[, j] == 1)
      temp$onset[disease_generated] <- 1
      temp$onset_age[disease_generated] <- disease_age[disease_generated, j]
      fit_model_cox(temp, cov_vector)
    })


    average_coefficients <- function(lst) {
      n <- 3
      lapply(seq_len(n), function(i) {
        components <- lapply(lst, `[[`, i)
        Reduce("+", components) / length(components)
      })
    }

    avg_coefficients <- average_coefficients(result_iteration)

    criteria1 <- norm(as.matrix((avg_coefficients[[1]] - old[[1]]) / old[[1]]))
    criteria2 <- norm(as.matrix((avg_coefficients[[2]] - old[[2]]) / old[[2]]))
    criteria3 <- norm(as.matrix((avg_coefficients[[2]] - old[[2]]) / old[[2]]))
    criteria <- criteria1+criteria2+criteria3


    lp.12 <- covar %*% avg_coefficients[[1]]
    lp.13 <- covar %*% avg_coefficients[[2]]
    lp.23 <- covar %*% avg_coefficients[[3]]

    print(avg_coefficients)
    old <- avg_coefficients

    k <- k+1
    message(sprintf(
      "Iteration %d — Total criteria: %.5f (c1: %.5f, c2: %.5f,  c3: %.5f)",
      k, criteria, criteria1, criteria2, criteria3
    ))
  }


  all_fits <-  parallel::mclapply(1:m, function(j) {
    temp <- scheme_data
    disease_generated <- which(disease_status[, j] == 1)
    temp$onset[disease_generated] <- 1
    temp$onset_age[disease_generated] <- disease_age[disease_generated, j]
    fit_model(temp, cov_vector, clock_assumption, distribution)
  })



  # Compute averaged parameters
  averaged_params <- averaging_params(all_fits)

  # Bootstrap

  bias_hat = NULL
  B <- 100
  m_boot <- 5

  if(boot == TRUE){
    boots <- bias_bootstrap(
      scheme_data = scheme_data,
      disease_status = disease_status,
      disease_age = disease_age,
      cov_vector = cov_vector,
      clock_assumption = clock_assumption,
      distribution = distribution,
      averaging_params = averaging_params,
      B = B, m_boot = m_boot, seed = 42
    )

    boot_mean <- Reduce("+", boots) / length(boots)
    bias_hat <- boot_mean - averaged_params
  }


  return(list(all_fits, averaged_params,bias_hat ))
}



fit_one_seed_mipd_iter <- function(seed,cov_vector,path,out_dir) {
  sim_path <- file.path(path,
                        sprintf("simulation_ready_%03d.rds", seed)
  )
  if (!file.exists(sim_path)) {
    message(sprintf("Missing file for seed %d -> skip", seed))
    return(NULL)
  }

  d <- readRDS(sim_path)
  panel_data <- d[[2]]

  m <- 20
  clock_assumption <- "forward"
  distribution <- "gompertz"
  inner_cores <- 1

  eps <- 1e-3
  max_iter <- 5


  res <- tryCatch(
    run_imputation_iterative(panel_data, m, cov_vector, clock_assumption, distribution, inner_cores, max_iter, eps),
    error = function(e) NULL
  )


  out_seed_dir <- file.path(out_dir)
  dir.create(out_seed_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(res[[1]],  file.path(out_seed_dir, sprintf("seed_%03d_model_mipd.rds",  seed)))


  list(params = res[[2]])
}


# mipd_iter bc

fit_one_seed_mipd_iter_bc <- function(seed,cov_vector,path,out_dir,boot = TRUE) {
  sim_path <- file.path(path,
                        sprintf("simulation_ready_%03d.rds", seed)
  )
  if (!file.exists(sim_path)) {
    message(sprintf("Missing file for seed %d -> skip", seed))
    return(NULL)
  }

  d <- readRDS(sim_path)
  panel_data <- d[[2]]

  m <- 20
  clock_assumption <- "forward"
  distribution <- "gompertz"
  inner_cores <- 20
  eps <- 1e-3
  max_iter <- 10


  res <- tryCatch(
    run_imputation_iterative(panel_data, m, cov_vector, clock_assumption, distribution, inner_cores, max_iter, eps, boot =TRUE),
    error = function(e) NULL
  )


  out_seed_dir <- file.path(out_dir)
  dir.create(out_seed_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(res[[1]],  file.path(out_seed_dir, sprintf("seed_%03d_model_mipd.rds",  seed)))
  saveRDS(res[[3]],  file.path(out_seed_dir, sprintf("seed_%03d_bias_mipd.rds",  seed)))


  list(params = res[[2]])
}

