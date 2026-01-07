simulate_panel_data <- function(n_pats, prevalence, max_sim_age, model_type, followup, nseed,rate,event_rate,shape_fixed,rate_fixed) {
  
  set.seed(123)
  
  # Step 1: Simulate disease progression
  simulation_data <- function(n_pats, model_type, prevalence, covariates, nseed) {
    define_coefficients <-  function(rate, covariates){
      coefs_01 <- list(
        shape = matrix(shape_fixed[1], nrow = 1),
        rate = matrix(c(event_rate[rate], covariates[1,]), nrow = 1)
      )
      colnames(coefs_01$rate) <- c("(Intercept)", "cov1", "cov2", "cov3")
      colnames(coefs_01$shape) <- c("(Intercept)")
      
      coefs_02 <- list(
        shape = matrix(shape_fixed[2], nrow = 1),
        rate = matrix(c(rate_fixed[1], covariates[2,]), nrow = 1)
      )
      colnames(coefs_02$rate) <- c("(Intercept)", "cov1", "cov2", "cov3")
      colnames(coefs_02$shape) <- c("(Intercept)")
      
      coefs_12 <- list(
        shape = matrix(shape_fixed[3], nrow = 1),
        rate = matrix(c(rate_fixed[2], covariates[3,]), nrow = 1)
      )
      colnames(coefs_12$rate) <- c("(Intercept)", "cov1", "cov2", "cov3")
      colnames(coefs_12$shape) <- c("(Intercept)")
      
      return(list(coefs_01,coefs_02,coefs_12))
    }
    covariates_distr <- function(n_obs, prevalence){
      cov1 <- rbinom(n=n_obs, size=1, prob=prevalence)
      cov2 <- rbinom(n=n_obs, size=1, prob=prevalence)
      cov3 <- rnorm(n=n_obs, mean =1 , sd =0.3 )
      x1 <- runif(n_obs, min = 60, max = 90)
      Intercept <- matrix(1,n_pats,1)
      cov_df <- data.frame(Intercept, x1, cov1, cov2, cov3)
      colnames(cov_df)[1]="(Intercept)"
      return(cov_df=cov_df)
    }
    create_input_data<- function(n_pats, prevalence,nseed){
      patients <- covariates_distr(n_pats*nseed,prevalence)
      patients$patient_id=1:(n_seeds*n_pats)
      patients$id=rep(1:n_pats,times=n_seeds)
      patients$grp_id=rep(1:nseed,each=n_pats)
      strategies <- data.frame(strategy_id = 1)
      hesim_obj <- hesim_data(strategies = strategies,
                              patients = patients)
      input_data <- hesim::expand(hesim_obj, by = c("strategies", "patients"), times = NULL )
      return(list(input_data[, ],patients))
    }
    
    dat <- create_input_data(n_pats, prevalence,nseed)
    input_data <-dat[[1]]
    entry_matrix <- matrix(0,n_pats,4)
    entry_matrix[,1] <- 1
    colnames(entry_matrix) <- c("x1","cov1", "cov2", "cov3")
    fic_tmat <- mstate::transMat(x = list(c(2),c(3, 4),c(4),c()), names = c("Zero","Disease-free","Disease", "Death"))
    det_model <-  params_surv( coefs = list(est = entry_matrix), dist="fixed")
    det_model$n_samples <- 1
    
    coef <- define_coefficients(rate, covariates)
    
    object <- params_surv_list(
      det_model,
      params_surv(coefs = coef[[1]], dist = "gompertz"),
      params_surv(coefs = coef[[2]], dist = "gompertz"),
      params_surv(coefs = coef[[3]], dist = "gompertz")
    )
    
    dismod <- create_IndivCtstmTrans(
      object,
      input_data,
      trans_mat = fic_tmat,
      clock = model_type,
      start_age = 0
    )
    return(list(dismod,dat[[2]]))
  }
  
  sim_data <- simulation_data(n_pats, model_type, prevalence, covariates, nseed)
  true_data <- sim_data[[1]]$sim_disease(max_t = max_sim_age, max_age = max_sim_age)
  true_data_original <- true_data
  true_data[, c("sample","strategy_id") := NULL]
  
  
  if (any(true_data$time_stop == max_sim_age & true_data$to < 4)){
    message("Solving death times")
    print(sum((true_data$time_stop == max_sim_age & true_data$to < 4)))
  }
  
  
  
  true_data <- true_data %>%
    filter(!(from == 1 & to == 2)) %>%
    mutate(
      to   = ifelse(to == 3, 2, 3),
      from = ifelse(from == 2, 1, 2)
    )
  
  if (any(true_data$time_start>true_data$time_stop))
    stop("Error with start-stop times")
  
  
  true_data <- true_data %>%
    left_join(sim_data[[2]]) %>%
    mutate(patient_id=id) %>%
    dplyr::select(-id) %>%
    rename(seed=grp_id)
  
  
  # Step 2: Right censoring
  apply_right_censoring <- function(temp,nseed){
    #right_censoring <- runif(n = n_pats, min = 0, max = 20)
    right_censoring <- rnorm(n = n_pats*nseed, mean=18, sd=3)
    rc_id <- tibble(right_censoring = right_censoring, seed=rep(1:nseed,each=n_pats),
                    patient_id=rep(1:n_pats,times=nseed))
    temp <- temp%>%
      left_join(rc_id)
    temp$right_censoring <- temp$right_censoring + temp$age
    temp$onset_age_c <- pmin(temp$onset_age, temp$death_time, temp$right_censoring, na.rm=T)
    temp$death_time_c <- pmin(temp$death_time, temp$right_censoring)
    temp$onset_c <- ifelse(temp$onset_age_c==temp$right_censoring | temp$onset_age_c== temp$death_time_c, 0, temp$onset)
    temp$dead <- ifelse(temp$death_time_c==temp$right_censoring, 0, 1)
    return(temp)
  }
  
  print(nrow(true_data))
  # Step 3: Introducing onset and death indicator/time
  true_data$age <- true_data$time_start
  true_data <- true_data %>%
    group_by(patient_id,seed) %>%
    mutate(age = if (n() == 2) first(time_start) else age,
           death_time = max(time_stop)) %>%
    ungroup() %>%
    mutate(onset = if_else(from == 1 & to == 2, 1, 0),
           onset_age = if_else(from==1 & to==2, time_stop, NA_integer_ ))
  print(nrow(true_data))
  true_data <- true_data %>%
    group_by(patient_id,seed) %>%
    slice(1) %>%
    ungroup()
  
  print(nrow(true_data))
  true_data <- apply_right_censoring(true_data,nseed)
  print(nrow(true_data))
  true_data <- true_data %>%
    select(-from, -to, -time_start, -time_stop, -final)
  
  
  
  onset_after_rc <- true_data$patient_id[true_data$onset_c==1]
  original_onset <- true_data$patient_id[true_data$onset==1]
  length(intersect(onset_after_rc,original_onset))
  print(nrow(true_data))
  
  # Step 4: Applying observation scheme and accordingly updating quantities
  observation_scheme <- function(data, row){
    data$row <- 1:nrow(data)
    which_onset <- which(data$onset_c == 1)
    patient_data <- data[data$row==row,]
    avg_spacing <- 3
    n_visits <- ceiling(followup / avg_spacing)
    obs_scheme <- runif(n_visits - 1, min = avg_spacing - 0.5, max = avg_spacing + 0.5)
    visit_times <- c(0, cumsum(obs_scheme)) + patient_data$age
    
    censor_index <- which(visit_times >= patient_data$death_time_c)[1]
    censor_index <- ifelse(is.na(censor_index), length(visit_times), max(censor_index, 2))
    visit_times <- visit_times[1:(censor_index - 1)]
    
    if (patient_data$row %in% which_onset) {
      onset_index <- which(visit_times >= patient_data$onset_age_c)[1]
      patient_data$onset_c <- ifelse(patient_data$onset_age_c > max(visit_times), 0, 1)
      patient_data$onset_age_c <- ifelse(patient_data$onset_c == 1, visit_times[onset_index], patient_data$death_time_c)
    } else {
      patient_data$onset_c <- 0
      patient_data$onset_age_c <- patient_data$death_time_c
    }
    
    n_repeats <- length(visit_times)
    patient_data <- patient_data[rep(1, n_repeats), ]
    patient_data$age <- visit_times
    patient_data$visits <- 1:n_repeats
    patient_data$onset_c <- ifelse(patient_data$age >= patient_data$onset_age_c, 1, 0)
    
    return(patient_data)
  }
  
  panel_data <- lapply(1:nrow(true_data), function(i) observation_scheme(true_data, i))
  panel_dataset <- do.call(rbind, panel_data)
  
  return(list(true_data = true_data, panel_dataset = panel_dataset))
}