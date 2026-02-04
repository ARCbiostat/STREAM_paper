

run_models <- function(rate,models,covnames,pathnamein="Simulation/simulation_results",pathnameout="Simulation/models_results",cores=1){

  plan(multisession,workers=cores)
  sim_root <- file.path(
    pathnamein,
    rate
  )
  
  print(sim_root)

  print(sim_root)
  gt_root  <- paste0(pathnamein,"/ground_truth_params")

  seeds <- 1:200



  out_root <- file.path(
    pathnameout,
    rate
  )

print(out_root)
#   # # =====================  survival model ==========================

  model_tag <- "cox"
if(model_tag%in%models){
  coef_list <- future_lapply(seeds, function(s) {
    cat(sprintf("\n[SEED %03d] Cox start\n", s))
    out <- tryCatch(fit_one_seed_cox(s,sim_root,out_root,model_tag,covariate_names=covnames),
                    error = function(e) { message(sprintf("[SEED %03d] ERROR: %s", s, e$message)); NULL })
    if (is.null(out)) cat(sprintf("[SEED %03d] Cox failed\n", s)) else cat(sprintf("[SEED %03d] Cox done\n", s))
    out
  })

  coef_list <- Filter(Negate(is.null), coef_list)
  if (length(coef_list) == 0) stop("No successful Cox fits; check paths/seeds/data.")

  all_terms <- sort(unique(unlist(lapply(coef_list, names))))

  coef_mat <- sapply(coef_list, function(v) {
    out <- rep(NA_real_, length(all_terms)); names(out) <- all_terms
    out[names(v)] <- unname(v)
    out
  })
  rownames(coef_mat) <- all_terms

  cox_params <- rowMeans(coef_mat, na.rm = TRUE)
  cox_params <- as.matrix(cox_params, byrow=T)
  colnames(cox_params) <- "cov"
  rownames(cox_params) <- c(1,2,3)

  out_base <- file.path(out_root,
                        model_tag)
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

  saveRDS(cox_params, file.path(out_base, "cox_params_avg.rds"))

  print("ok Cox")
}

 # ==========================================================

 # ========================== msm ============================

  model_tag  <- "msm"
  if(model_tag%in%models){
    out_base <- file.path(out_root,
                          model_tag)
    dir.create(out_base, recursive = TRUE, showWarnings = FALSE)


    fits_list <-lapply(seeds, function(s) fit_one_seed_msm(s, age = FALSE,sim_root,out_base,covariate_names=covnames))
    fits_list <- Filter(function(x) {   !is.null(x$params) && all(!is.na(x$params)) }, fits_list)
    if (length(fits_list) == 0) stop("No simulations loaded; check paths/seeds.")


    params_only <- lapply(fits_list, function(x) x$params)

    gc()
    # ==========================================================


    print("ok msm")
  }


# # ========================== msm + age  ============================

  model_tag  <- "msm_age"

  if(model_tag %in% models){
    out_base <- file.path(out_root,
                          model_tag)
    dir.create(out_base, recursive = TRUE, showWarnings = FALSE)


    fits_list <- future_lapply(seeds, function(s) fit_one_seed_msm(s, age = TRUE,sim_root,out_base,covariate_names=covnames))
    fits_list <- Filter(function(x) {   !is.null(x$params) && all(!is.na(x$params)) }, fits_list)
    if (length(fits_list) == 0) stop("No simulations loaded; check paths/seeds.")

    gc()
    print("ok msm age")
  }


# # =========================================================
#
#
# # =========== flexsurv over real and panel data =============
  #

  model_tag  <- "flexsurv"

  if(model_tag%in%models){
    out_base <- file.path(out_root,
                          model_tag)
    dir.create(out_base, recursive = TRUE, showWarnings = FALSE)


    res_list <-future_lapply(seeds, function(x)fit_one_seed_gomp(x,sim_root,out_root,model_tag,covariate_names=covnames))
    res_list <- Filter(Negate(is.null), res_list)

    if (length(res_list) == 0) stop("No simulations loaded; check paths or seeds.")


    n_trans <- nrow(res_list[[1]]$real)
    param_names <- colnames(res_list[[1]]$real)


    avg_mat <- function(which_mat = c("real","panel")) {
      which_mat <- match.arg(which_mat)
      out <- matrix(NA_real_, nrow = n_trans, ncol = length(param_names),
                    dimnames = list(rownames(res_list[[1]]$real), param_names))
      for (i in 1:n_trans) {
        for (j in seq_along(param_names)) {
          vals <- sapply(res_list, function(x) x[[which_mat]][i, j])
          out[i, j] <- mean(vals, na.rm = TRUE)
        }
      }
      out
    }

    params_real  <- avg_mat("real")
    params_panel <- avg_mat("panel")

    saveRDS(params_real, file.path(out_base, "params_real_avg.rds"))
    saveRDS(params_panel, file.path(out_base, "params_panel_avg.rds"))

    gc()

    print("ok flex")
  }


# ==========================================================



# =========================  nhm  ===========================

  model_tag  <- "nhm"
  if(model_tag%in%models){
    out_base <- file.path(out_root,
                          model_tag)
    dir.create(out_base, recursive = TRUE, showWarnings = FALSE)


    fits_list <-future_lapply(seeds, function(s) fit_one_seed_nhm(s,sim_root,sim_root,out_base,covariate_names=covnames)
                         )
    fits_list <- Filter(function(x) {   !is.null(x$params) && all(!is.na(x$params)) }, fits_list)
    fits_list <- Filter(function(x) !isTRUE(x$model$singular), fits_list)


    if (length(fits_list) == 0) stop("No simulations loaded; check paths/seeds.")


    gc()

    # ==========================================================

    print("ok nhm")
  }

  # =========================  mipd_iterative  ===========================
#
model_tag  <- "mipd_iter"

  if(model_tag%in%models){
    out_base <- file.path(out_root,
                          model_tag)
    dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

    cov_vector<-"cov"

    fits_list <-future_lapply(seeds, function(s) fit_one_seed_mipd_iter(s, cov_vector, sim_root,out_base,covariate_names=covnames))
    fits_list <- Filter(function(x) {   !is.null(x$params) && all(!is.na(x$params)) }, fits_list)
    if (length(fits_list) == 0) stop("No simulations loaded; check paths/seeds.")


    params_only <- lapply(fits_list, function(x) x$params)


print("ok mipd iter")
  }


#=========================streams=====================================


model_tag  <- "streams"

if(model_tag%in%models){
  out_base <- file.path(out_root,
                        model_tag)
  dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
  
  cov_vector <- list(
    "0->1" = covnames,
    "0->2" = covnames,
    "1->2" = covnames
  )
  
  fits_list <-future_lapply(seeds, function(s) fit_one_seed_streams(s, cov_vector, sim_root,out_base,covariate_names=covnames))
  fits_list <- Filter(function(x) {   !is.null(x$params) && all(!is.na(x$params)) }, fits_list)
  if (length(fits_list) == 0) stop("No simulations loaded; check paths/seeds.")
}

print("I AM DONE.")

stopCluster()
gc()
}


