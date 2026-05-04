library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tidyr)
library(STREAMS)

source("fun/utils_stat.R")

# ===================== CONFIG =====================

term_rate  <- c("low")
seeds      <- 1:200
covnames <- c("cov1","cov2","cov3")

gt_root  <- "Simulation/simulation_results/ground_truth_params"

sim_root <- file.path(
  "Simulation/models_results_unif"
)
out_root <- file.path(
  "Simulation/model_estimates_unif"
)

model_tags <- c("flexsurv","msm_age","nhm","mipd_iter","streams","streams2")
# LOW

model_tags <- "msm_age"
rate <- term_rate[1]

# ----------------- main -----------------
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)



  ground_truth_base <- file.path(gt_root, rate,"params_scheme.rds")
  ground_truth <- readRDS(ground_truth_base)

  scheme_out_dir <- file.path(out_root, rate)
  dir.create(scheme_out_dir, recursive = TRUE, showWarnings = FALSE)

  for (tag in model_tags) {
    cat(" • Model tag:", tag, "\n")

    files_df <- list_model_files(sim_root, rate, tag)
    if (nrow(files_df) == 0) {
      warning(sprintf(" %s] no files found", tag))
      next
    }

    if (!is.null(seeds)) {
      files_df <- dplyr::filter(files_df, seed %in% seeds)
      if (nrow(files_df) == 0) {
        warning(sprintf(" %s] no matching seeds", tag))
        next
      }
    }

    per_seed <- purrr::pmap_dfr(
      files_df,
      function(filepath, seed, model_sub = NA_character_, kind =NULL) {
        obj <- readRDS(filepath)


        dispatch_tag <- if (tag == "flexsurv") "flexsurv" else tag
        covnames_all <- if (tag=="msm_age")  c(covnames,"age") else covnames
        ci <- try(
           safe_ci_call(dispatch_tag, obj,covnames_all),
           silent = TRUE)
        if (inherits(ci, "try-error")) {
          warning(sprintf("[  %s] seed %s: CI extraction failed",
                          tag, seed))
          return(tibble::tibble())
        }


        reported_model <- if (!is.na(model_sub)) model_sub else tag

        dplyr::mutate(ci, model = reported_model, seed = seed,rate=rate)
      }
    )

    if (nrow(per_seed) == 0) {
      warning(sprintf("[ %s] nothing extracted", tag))
      next
    }

    ic_summary <- summarize_across_seeds(per_seed)
    out_base <- file.path(out_root, rate, tag)
    dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
    saveRDS(ic_summary, file = file.path(out_base, "ic_summary.rds"))
  }




  # MEDIUM
  
  rate <- term_rate[2]
  
  # ----------------- main -----------------
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  
  
  
  ground_truth_base <- file.path(gt_root, rate,"params_scheme.rds")
  ground_truth <- readRDS(ground_truth_base)
  
  scheme_out_dir <- file.path(out_root, rate)
  dir.create(scheme_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (tag in model_tags) {
    cat(" • Model tag:", tag, "\n")
    
    files_df <- list_model_files(sim_root, rate, tag)
    if (nrow(files_df) == 0) {
      warning(sprintf(" %s] no files found", tag))
      next
    }
    
    if (!is.null(seeds)) {
      files_df <- dplyr::filter(files_df, seed %in% seeds)
      if (nrow(files_df) == 0) {
        warning(sprintf(" %s] no matching seeds", tag))
        next
      }
    }
    
    per_seed <- purrr::pmap_dfr(
      files_df,
      function(filepath, seed, model_sub = NA_character_, kind =NULL) {
        obj <- readRDS(filepath)
        
        
        dispatch_tag <- if (tag == "flexsurv") "flexsurv" else tag
        covnames_all <- if (tag=="msm_age")  c("age",covnames) else covnames
        ci <- try(
          safe_ci_call(dispatch_tag, obj,covnames_all),
          silent = TRUE)
        if (inherits(ci, "try-error")) {
          warning(sprintf("[  %s] seed %s: CI extraction failed",
                          tag, seed))
          return(tibble::tibble())
        }
        
        
        reported_model <- if (!is.na(model_sub)) model_sub else tag
        
        dplyr::mutate(ci, model = reported_model, seed = seed,rate=rate)
      }
    )
    
    if (nrow(per_seed) == 0) {
      warning(sprintf("[ %s] nothing extracted", tag))
      next
    }
    
    ic_summary <- summarize_across_seeds(per_seed)
    out_base <- file.path(out_root, rate, tag)
    dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
    saveRDS(ic_summary, file = file.path(out_base, "ic_summary.rds"))
  }
  

  # HIGH
  
  rate <- term_rate[3]
  
  # ----------------- main -----------------
  dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
  
  
  
  ground_truth_base <- file.path(gt_root, rate,"params_scheme.rds")
  ground_truth <- readRDS(ground_truth_base)
  
  scheme_out_dir <- file.path(out_root, rate)
  dir.create(scheme_out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (tag in model_tags) {
    cat(" • Model tag:", tag, "\n")
    
    files_df <- list_model_files(sim_root, rate, tag)
    if (nrow(files_df) == 0) {
      warning(sprintf(" %s] no files found", tag))
      next
    }
    
    if (!is.null(seeds)) {
      files_df <- dplyr::filter(files_df, seed %in% seeds)
      if (nrow(files_df) == 0) {
        warning(sprintf(" %s] no matching seeds", tag))
        next
      }
    }
    
    per_seed <- purrr::pmap_dfr(
      files_df,
      function(filepath, seed, model_sub = NA_character_, kind =NULL) {
        obj <- readRDS(filepath)
        
        
        dispatch_tag <- if (tag == "flexsurv") "flexsurv" else tag
        covnames_all <- if (tag=="msm_age")  c("age",covnames) else covnames
        ci <- try(
          safe_ci_call(dispatch_tag, obj,covnames_all),
          silent = TRUE)
        if (inherits(ci, "try-error")) {
          warning(sprintf("[  %s] seed %s: CI extraction failed",
                          tag, seed))
          return(tibble::tibble())
        }
        
        
        reported_model <- if (!is.na(model_sub)) model_sub else tag
        
        dplyr::mutate(ci, model = reported_model, seed = seed,rate=rate)
      }
    )
    
    if (nrow(per_seed) == 0) {
      warning(sprintf("[ %s] nothing extracted", tag))
      next
    }
    
    ic_summary <- summarize_across_seeds(per_seed)
    out_base <- file.path(out_root, rate, tag)
    dir.create(out_base, showWarnings = FALSE, recursive = TRUE)
    saveRDS(ic_summary, file = file.path(out_base, "ic_summary.rds"))
  }
  

