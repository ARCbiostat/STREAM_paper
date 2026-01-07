

run_models <- function(size,rate,scheme_id,models,pathnamein="Simulation/simulation_results",pathnameout="Simulation/models_results"){

  sim_root <- file.path(
    pathnamein,
    size,
    paste0("cov_scheme_",scheme_id),
    rate
  )

  print(sim_root)
  gt_root  <- paste0(pathnamein,"/ground_truth_params")

  seeds <- 1:200



  out_root <- file.path(
    pathnameout,
    size,
    paste0("cov_scheme_",scheme_id),
    rate
  )


#   # # =====================  survival model ==========================

  model_tag <- "cox"
if(model_tag%in%models){
  coef_list <- mclapply(seeds, function(s) {
    cat(sprintf("\n[SEED %03d] Cox start\n", s))
    out <- tryCatch(fit_one_seed_cox(s,sim_root,out_root,model_tag),
                    error = function(e) { message(sprintf("[SEED %03d] ERROR: %s", s, e$message)); NULL })
    if (is.null(out)) cat(sprintf("[SEED %03d] Cox failed\n", s)) else cat(sprintf("[SEED %03d] Cox done\n", s))
    out
  }
  ,  mc.cores = 10)

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


    fits_list <-mclapply(seeds, function(s) fit_one_seed_msm(s, age = FALSE,sim_root,out_base), mc.cores = 10)
    fits_list <- Filter(function(x) {   !is.null(x$params) && all(!is.na(x$params)) }, fits_list)
    if (length(fits_list) == 0) stop("No simulations loaded; check paths/seeds.")


    params_only <- lapply(fits_list, function(x) x$params)

    row_names <- rownames(params_only[[1]])
    col_names <- colnames(params_only[[1]])

    arr <- simplify2array(lapply(params_only, function(M) {
      M <- as.matrix(M)
      dimnames(M) <- list(row_names, col_names)
      M
    }))

    msm_params <- apply(arr, c(1, 2), function(x) mean(x, na.rm = TRUE))

    saveRDS(msm_params, file.path(out_base, "msm_params_avg.rds"))

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


    fits_list <- mclapply(seeds, function(s) fit_one_seed_msm(s, age = TRUE,sim_root,out_base), mc.cores = 10)
    fits_list <- Filter(function(x) {   !is.null(x$params) && all(!is.na(x$params)) }, fits_list)
    if (length(fits_list) == 0) stop("No simulations loaded; check paths/seeds.")


    # params_only <- lapply(fits_list, function(x) x$params)
    #
    # row_names <- rownames(params_only[[1]])
    # col_names <- colnames(params_only[[1]])
    #
    # arr <- simplify2array(lapply(params_only, function(M) {
    #   M <- as.matrix(M)
    #   dimnames(M) <- list(row_names, col_names)
    #   M
    # }))

    # msm_age_params <- apply(arr, c(1, 2), function(x) mean(x, na.rm = TRUE))
    # msm_age_params <- round(msm_age_params,6)
    # saveRDS(msm_age_params, file.path(out_base, "msm_age_params_avg.rds"))

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


    res_list <-lapply(seeds, function(x)fit_one_seed_gomp(x,sim_root,out_root,model_tag))
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


    fits_list <-mclapply(seeds, function(s) fit_one_seed_nhm(s,sim_root,sim_root,out_base)
                         ,  mc.cores = 10)
    fits_list <- Filter(function(x) {   !is.null(x$params) && all(!is.na(x$params)) }, fits_list)
    fits_list <- Filter(function(x) !isTRUE(x$model$singular), fits_list)


    if (length(fits_list) == 0) stop("No simulations loaded; check paths/seeds.")


    params_only <- lapply(fits_list, function(x) x$params)

    row_names <- rownames(params_only[[1]])
    col_names <- colnames(params_only[[1]])

    arr <- simplify2array(lapply(params_only, function(M) {
      M <- as.matrix(M)
      dimnames(M) <- list(row_names, col_names)
      M
    }))

    # nhm_params <- apply(arr, c(1, 2), function(x) mean(x, na.rm = TRUE))
    #
    # saveRDS(nhm_params, file.path(out_base, "nhm_params_avg.rds"))


    gc()

    # ==========================================================

    print("ok nhm")
  }


# =========================  mipd  ===========================


# model_tag  <- "mipd"
#
# out_base <- file.path(out_root,
#                       model_tag)
# dir.create(out_base, recursive = TRUE, showWarnings = FALSE)
#
# cov_vector<-"cov"
# fits_list <-mclapply(seeds, function(s) fit_one_seed_mipd(s, cov_vector,sim_root,out_base)
#                      ,  mc.cores = 4)
# fits_list <- Filter(function(x) {   !is.null(x$params) && all(!is.na(x$params)) }, fits_list)
# if (length(fits_list) == 0) stop("No simulations loaded; check paths/seeds.")
#
#
# params_only <- lapply(fits_list, function(x) x$params)
#
# row_names <- rownames(params_only[[1]])
# col_names <- colnames(params_only[[1]])
#
# arr <- simplify2array(lapply(params_only, function(M) {
#   M <- as.matrix(M)
#   dimnames(M) <- list(row_names, col_names)
#   M
# }))
#
# mipd_params <- apply(arr, c(1, 2), function(x) mean(x, na.rm = TRUE))
#
# saveRDS(mipd_params, file.path(out_base, "mipd_params_avg.rds"))




  #==============================================================





  # =========================  mipd_iterative  ===========================
#
model_tag  <- "mipd_iter"

  if(model_tag%in%models){
    out_base <- file.path(out_root,
                          model_tag)
    dir.create(out_base, recursive = TRUE, showWarnings = FALSE)

    cov_vector<-"cov"

    fits_list <-mclapply(seeds, function(s) fit_one_seed_mipd_iter(s, cov_vector, sim_root,out_base),  mc.cores = 1)
    fits_list <- Filter(function(x) {   !is.null(x$params) && all(!is.na(x$params)) }, fits_list)
    if (length(fits_list) == 0) stop("No simulations loaded; check paths/seeds.")


    params_only <- lapply(fits_list, function(x) x$params)

    row_names <- rownames(params_only[[1]])
    col_names <- colnames(params_only[[1]])

    arr <- simplify2array(lapply(params_only, function(M) {
      M <- as.matrix(M)
      dimnames(M) <- list(row_names, col_names)
      M
    }))

    # mipd_iter_params <- apply(arr, c(1, 2), function(x) mean(x, na.rm = TRUE))
    #
    # saveRDS(mipd_iter_params, file.path(out_base, "mipd_iter_params_avg.rds"))
print("ok mipd iter")
  }


#==============================================================



#================== Teacher-student traiing ===========================

# model_tag  <- "teach_stud"
# cov_vector <- c ("cov")
# cond_vector   <- c("cov", "age", "interval" )
# inner_cores <- 10
#
# # mclapply(seeds, function(s) fit_teach_stud(seed, cov_vector, cond_vector, sim_root, inner_cores, model_tag),  mc.cores = 10)
#
# for (seed in seeds){
#   fit_teach_stud(seed, cov_vector, cond_vector, sim_root, inner_cores, model_tags, scheme_id, rate)
#   print(seed)
# }
#


gc()


print("I AM DONE.")


}


run_model_cvae <- function(size,rate,scheme_id,pathnamein="Simulation/simulation_results",pathnameout="Simulation/models_results",boot){

  sim_root <- file.path(
    pathnamein,
    size,
    paste0("cov_scheme_",scheme_id),
    rate
  )

  print(sim_root)
  gt_root  <- paste0(pathnamein,"/ground_truth_params")

  seeds <- 1:200

  model_tag  <- "teach_stud"

  out_root <- file.path(
    pathnameout,
    size,
    paste0("cov_scheme_",scheme_id),
    rate,
    model_tag
  )

  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

  cov_vector <- c ("cov")
  cond_vector   <- c("cov", "age", "interval" )
  inner_cores <- 1

  # mclapply(seeds, function(s) fit_teach_stud(seed, cov_vector, cond_vector, sim_root, inner_cores, model_tag),  mc.cores = 10)

  for (seed in seeds){
    fit_teach_stud(seed, cov_vector, cond_vector, sim_root, inner_cores, model_tags, scheme_id, rate, pathout=out_root,boot)
    print(seed)
  }


}


