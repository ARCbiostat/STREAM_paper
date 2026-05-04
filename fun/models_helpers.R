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

fit_one_seed_cox <- function(seed,path,out_dir,model_tag,covariate_names) {
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

  fit <- fit_cox(temp,covariate_names=covariate_names)

  coefs <- stats::coef(fit)


  out_dir <- file.path(out_dir,
                       model_tag)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


  saveRDS(fit,   file.path(out_dir, sprintf("seed_%03d_model_cox.rds", seed)))


  coefs
}

# flexsurv

prep_long <- function(df, tmat,covariate_names) {
  dl <- mstate::msprep(
    data   = df,
    trans  = tmat,
    time   = c(NA, "onset_age", "death_time"),
    status = c(NA, "onset", "dead"),
    keep   = c("age",covariate_names),
    id     = "patient_id"
  )
  dl$Tstart[dl$trans < 3] <- dl$Tstart[dl$trans < 3] + dl$age[dl$trans < 3]
  dl %>% mutate(time = Tstop - Tstart)
}

fit_one_seed_gomp <- function(seed,path,out_dir,model_tag,covariate_names) {
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

  data_long_real  <- prep_long(real_data,  tmat,covariate_names)
  data_long_panel <- prep_long(panel_collapsed, tmat,covariate_names)

  n_trans <- max(tmat, na.rm = TRUE)


  fits_real  <- vector("list", n_trans)
  fits_panel <- vector("list", n_trans)
  for (i in 1:n_trans) {
    sub_r <- subset(data_long_real,  trans == i)
    sub_p <- subset(data_long_panel, trans == i)

    fits_real[[i]] <- if (nrow(sub_r) == 0 || sum(sub_r$status) == 0) NA else
      tryCatch(flexsurvreg(as.formula(paste("Surv(Tstart, Tstop, status) ~", paste(covariate_names, collapse = " + "))), data = sub_r, dist = "gompertz"),
               error = function(e) NA)

    fits_panel[[i]] <- if (nrow(sub_p) == 0 || sum(sub_p$status) == 0) NA else
      tryCatch(flexsurvreg(as.formula(paste("Surv(Tstart, Tstop, status) ~", paste(covariate_names, collapse = " + "))), data = sub_p, dist = "gompertz"),
               error = function(e) NA)
  }


  get_param_names <- function(fits,covariate_names) {
    for (k in seq_along(fits)) {
      if (is.list(fits[[k]]) && !is.na(fits[[k]])[1]) return(names(fits[[k]]$coefficients))
    }

    return(c("shape","rate",covariate_names))
  }
  param_names <- get_param_names(fits_panel,covariate_names)


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
  panel_data_death <- df%>% filter(dead==1 ) %>% group_by(patient_id) %>% mutate(age=death_time) %>% distinct(patient_id,.keep_all = T)
  df2 <- rbind(df,panel_data_death) %>% group_by(patient_id) %>%  mutate(state=case_when(max(onset)==1 & onset_age<=age& age!=death_time~2 ,
                                                              death_time==age & dead==1 ~3,
        TRUE~1)) %>% arrange(patient_id)
  return(df2)
}

fit_msm <- function(df,covariate_names) {
  Q <- rbind(c(0, 1, 1),
             c(0, 0, 1),
             c(0, 0, 0))

  fit <- tryCatch(
    msm(state ~ age,
        subject    = patient_id,
        data       = df,
        qmatrix    = Q,
        covariates = as.formula(paste("~",paste(covariate_names, collapse = " + "))),
        censor     = 99,
        control    = list(fnscale = 1000, maxit = 1000),
        deathexact = TRUE),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 1+length(covariate_names),
                  dimnames = list(c("1","2","3"),
                                  c("rate",covariate_names)))
    return(list(model = NULL, params = mat))
  }

  v <- fit$estimates
  if (length(v) < 6) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 1+length(covariate_names),  byrow = FALSE,
                  dimnames = list(c("1","2","3"),
                                  c("rate",covariate_names)))
  } else {
    mat <- matrix(v, nrow = 3, ncol = 1+length(covariate_names), byrow = FALSE,
                  dimnames = list(c("1","2","3"),
                                  c("rate",covariate_names)))
  }

  list(model = fit, params = mat)
}

fit_msm_age <- function(df,seed,out_dir,covariate_names) {

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
        covariates = as.formula(paste("~",paste(covariate_names, collapse = " + "),"+age")) ,
        control    = list(fnscale = 1000, maxit = 1000),
        censor = 99,
        deathexact = TRUE),
    error = function(e) NULL
  )

  if (is.null(fit)) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 3,
                  dimnames = list(c("1","2","3"),
                                  c("rate",covariate_names, "age")))
    return(list(model = NULL, params = mat))
  }

  v <- fit$estimates
  if (length(v) < 9) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 3, byrow = F,
                  dimnames = list(c("1","2","3"),
                                  c("rate",covariate_names, "age")))
  } else {
    mat <- matrix(v, nrow = 3,ncol=2+length(covariate_names), byrow = F,
                  dimnames = list(c("1","2","3"),
                                  c("rate",covariate_names, "age")))
  }

  list(model = fit, params = mat)
}

fit_one_seed_msm <- function(seed, age= FALSE,path,out_dir,covariate_names) {
  sim_path <- file.path(path,
                        sprintf("simulation_ready_%03d.rds", seed)
  )
  if (!file.exists(sim_path)) {
    message(sprintf("Missing file for seed %d -> skip", seed))
    return(NULL)
  }

  d <- readRDS(sim_path)
  panel_data <- d[[2]]
  temp_panel <- prepare_msm(panel_data)
  out_base <- out_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (age == FALSE){
    res <- fit_msm(temp_panel,covariate_names)
    saveRDS(res$model,  file.path(out_base, sprintf("seed_%03d_model_msm.rds",  seed)))
  }else{
    res <- fit_msm_age(temp_panel, seed,out_base,covariate_names)
    saveRDS(res$model,  file.path(out_base, sprintf("seed_%03d_model_msm_age.rds",  seed)))
  }

  list( params = res$params)
}

# nhm

fit_nhm <- function(df,covariate_names) {
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

  covm <- rep(list(tmat_1), length(covariate_names))
  covm <- lapply(seq_len(length(covariate_names)), function(i) (covm[[i]] + (i-1)*3)*matrix(c(0,1,1,
                                                                                              0,0,1,
                                                                                              0,0,0),byrow=T,ncol=3,nrow=3))
  names(covm) <- covariate_names
  fit <- tryCatch({
    object_nhm <- model.nhm(
      state ~ age,
      subject       = patient_id,
      data          = df,
      trans         = tmat_1,
      nonh          = tmat_1,
      type          = "gompertz",
      covariates    = covariate_names,
      covm          = covm,
      death         = TRUE,
      death.states  = c(3)
    )

    model_nhm <- nhm(
      object_nhm,
      gen_inits   = TRUE,
      score_test  = FALSE,
      control     = nhm.control(splits = split_points, obsinfo = FALSE)
    )

    model_nhm
  },
  error = function(e) NULL
  )


  if (is.null(fit)) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 2+length(covariate_names),
                  dimnames = list(c("1","2","3"),
                                  c("rate", "shape", covariate_names)))
    return(list(model = NULL, params = mat))
  }

  v <- fit$par
  if (length(v) < 9) {
    mat <- matrix(NA_real_, nrow = 3, ncol = 2+length(covariate_names),
                  dimnames = list(c("1","2","3"),
                                  c("rate", "shape", covariate_names)))
  } else {
    mat <- matrix(v[1:9], nrow = 3, ncol = 2+length(covariate_names),
                  dimnames = list(c("1","2","3"),
                                  c("rate", "shape", covariate_names)))
  }

  list(model = fit, params = mat)
}

fit_one_seed_nhm <- function(seed, age= FALSE,path,out_dir,covariate_names) {
  sim_path <- file.path(path,
                        sprintf("simulation_ready_%03d.rds", seed))
  if (!file.exists(sim_path)) {
    message(sprintf("Missing file for seed %d -> skip", seed))
    return(NULL)
  }

  d <- readRDS(sim_path)
  panel_data <- d[[2]]
  temp_panel <- prepare_msm(panel_data)

  res <- fit_nhm(temp_panel,covariate_names)
  out_base <- out_dir
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  saveRDS(res$model,  file.path(out_base, sprintf("seed_%03d_model_nhm.rds",  seed)))


  list( params = res$params, model = res$model)
}

# mipd

fit_one_seed_mipd_iter <- function(seed,cov_vector,path,out_dir,covariate_names) {
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
  panel_data %<>% select(patient_id,
                         seed,
                         age,
                         dead,
                         death_time,
                         onset,
                         onset_age,
                         visits,
                         any_of(covariate_names))
  
  


res <- tryCatch(
  run_mipd(panel_data,
m = m,
cov_vector = covariate_names,
clock_assumption = clock_assumption,
distribution = distribution,
inner_cores = 1,
max_iter = max_iter,
eps = eps),
    error = function(e) NULL
  )


  out_seed_dir <- file.path(out_dir)
  dir.create(out_seed_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(res,  file.path(out_seed_dir, sprintf("seed_%03d_model_mipd.rds",  seed)))


  list(params = res)
}



fit_one_seed_streams <- function(seed,cov_vector,path,out_dir,covariate_names) {
  sim_path <- file.path(path,
                        sprintf("simulation_ready_%03d.rds", seed)
  )
  if (!file.exists(sim_path)) {
    message(sprintf("Missing file for seed %d -> skip", seed))
    return(NULL)
  }
  
  d <- readRDS(sim_path)
  panel_data <- d[[2]]
  
  
  panel_data %<>% select(patient_id,
                         seed,
                         age,
                         dead,
                         death_time,
                         onset,
                         onset_age,
                         visits,
                         any_of(covariate_names))
  
  
  
  
  res <- #tryCatch(
    est <- run_streams(
      data = panel_data,
      cov_vector = cov_vector,
      python = Sys.which("python"),
      pu_args = list(verbose = TRUE),
      m=20
    #),
   # error = function(e) NULL
  )
  
  
  out_seed_dir <- file.path(out_dir)
  dir.create(out_seed_dir, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(res,  file.path(out_seed_dir, sprintf("seed_%03d_model_streams.rds",  seed)))
  
  
  list(params = res)
}


fit_one_seed_streams2 <- function(seed,cov_vector,path,out_dir,covariate_names) {
  sim_path <- file.path(path,
                        sprintf("simulation_ready_%03d.rds", seed)
  )
  if (!file.exists(sim_path)) {
    message(sprintf("Missing file for seed %d -> skip", seed))
    return(NULL)
  }
  
  d <- readRDS(sim_path)
  panel_data <- d[[2]]
  
  
  panel_data %<>% select(patient_id,
                         seed,
                         age,
                         dead,
                         death_time,
                         onset,
                         onset_age,
                         visits,
                         any_of(covariate_names))
  
  
  
  
  res <- #tryCatch(
    est <- run_streams2(
      data = panel_data,
      cov_vector = cov_vector,
      python = Sys.which("python"),
      pu_args = list(verbose = TRUE),
      m=20
      #),
      # error = function(e) NULL
    )
  
  
  out_seed_dir <- file.path(out_dir)
  dir.create(out_seed_dir, recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(res,  file.path(out_seed_dir, sprintf("seed_%03d_model_streams.rds",  seed)))
  
  
  list(params = res)
}

