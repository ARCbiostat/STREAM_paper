library(dplyr)
library(haven)
library(mstate)
library(flexsurv)
library(lubridate)
library(hesim)
library(here)
library(data.table)
library(tidyr)
library(tibble)
library(readr)


source("fun/simulate_panel_data.R")

base_dir <- "Simulation/simulation_results"
dir.create(base_dir, showWarnings = FALSE)

event_rate_labels <- c("low", "medium", "high")


base_dir <- "Simulation/simulation_results"
dir.create(base_dir, showWarnings = FALSE)

n_seeds <- 200
n_pats <- 2000
prevalence <- 0.5
model_type <- "forward"
followup <- 20
covariates <- t(rbind(
  c( 0.0, 0.4, 0.4),
  c( 0.7, 0.0, 0.0),
  c( 0.2, 0.2, 0.5)
))
rownames(covariates) <- c("B01","B02","B12")
colnames(covariates) <- c ("cov1", "cov2", "cov3")

rate <- c(-15, -14.15, -14)
event_rate_labels <- c("low", "medium", "high")
names(rate) <- event_rate_labels
shape_fixed     <- c("0->1" = 0.13,"0->2" = 0.11,  "1->2" = 0.075)
rate_fixed <- c( "0->2" = -12.5, "1->2" = -8.0)
transitions     <- c("0->1","0->2","1->2")

output_root <- "Simulation/simulation_results/ground_truth_params"

max_sim_age = 125

subfolder <- file.path(base_dir)

for (term in event_rate_labels) {
  outdir <- file.path(output_root, term)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  params <- cbind(
    shape = c(shape_fixed["0->1"], shape_fixed["0->2"], shape_fixed["1->2"]),
    rate  =c(rate[term], rate_fixed["0->2"], rate_fixed["1->2"]),
    cov1   = unname(covariates[,1]),
    cov2   = unname(covariates[,2]),
    cov3   = unname(covariates[,3])
    
  )
  rownames(params) <- transitions
    
    saveRDS(params, file.path(outdir, "params_scheme.rds"))
  }

cat(sprintf("Running simulation"))
output_root <- "Simulation/simulation_results"
for (term in event_rate_labels) {
  outdir <- file.path(output_root, term)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
res <-
  simulate_panel_data(
    n_pats      = n_pats,
    prevalence  = prevalence,
    max_sim_age = max_sim_age,
    model_type  = model_type,
    followup    = followup,
    nseed        = n_seeds,
    rate        =  term,
    event_rate  = rate,
    rate_fixed=rate_fixed,
    shape_fixed=shape_fixed
  )

ok <- is.list(res) &&
  all(c("true_data","panel_dataset") %in% names(res)) &&
  nrow(res$panel_dataset) > 0
if (!ok) { message("SKIP (empty or malformed result)"); next }

# save
result_file <- file.path(outdir, "simulation_data.rds")
saveRDS(res, result_file)
cat("saved ✓\n")

clean_cols <- function(df) {
  df %>%
    mutate(
      onset_real = onset,
      onset_age_real = onset_age,
      death_time_real = death_time,
      onset      = onset_c,
      onset_age  = onset_age_c,
      death_time = death_time_c
    ) %>%
    select(-onset_c, -onset_age_c, -death_time_c, -right_censoring)
}


get_subset_ids <- function(data, size) {
  n1 <- round(size * mean(data$onset == 1))
  n0 <- size - n1
  onset_1_ids <- sample(data$patient_id[data$onset == 1], n1)
  onset_0_ids <- sample(data$patient_id[data$onset == 0], n0)
  c(onset_1_ids, onset_0_ids)
}

lapply(1:n_seeds,function(s){
  print(result_file)
  
  if (!file.exists(result_file)) {
    message(cat("Missing file-> skip"))
    next
  }
  
  data <- readRDS(result_file)
  real_data  <- clean_cols(data[[1]]) %>% filter(seed==s)
  panel_data <- clean_cols(data[[2]]) %>% filter(seed==s)
  
  
  
  data_out <- list(real_data, panel_data)
  
  result_file <- file.path(outdir, sprintf("simulation_ready_%03d.rds", s))
  saveRDS(data_out, result_file)
  print(sprintf("Saved simulation: event_rate=%s, seed=%d",
                 term,s))
})


}








