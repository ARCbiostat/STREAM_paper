library(msm)
library(survival)
library(flexsurv)
library(dplyr)
library(data.table)
library(nhm)
library(here)
library(devtools)
library(truncnorm)
library(future)
library(future.apply)
library(Rcpp)
library(MIPD)
library(STREAMS)
library(magrittr)

source("fun/models_helpers.R")
source("fun/running_models_fun.R")

rates <- c("low","medium","high")
cores <- 4

# LOW EVENT RATES
system.time(run_models(rates[1],models=c("flexsurv","msm","msm_age","nhm","mipd_iter","streams"),covnames = c("cov1","cov2","cov3"),cores=cores))

print("finished low rates")
# MEDIUM EVENT RATES
gc()
system.time(run_models(rates[2],models=c("flexsurv","msm_age","nhm","mipd_iter","streams"),covnames = c("cov1","cov2","cov3"),cores=cores))
print("finished med rates")

# HIGH EVENT RATES
gc()
system.time(run_models(rates[3],scheme_id = 1,models=c("flexsurv","msm","msm_age","nhm","mipd_iter"),covnames = c("cov"),cores=cores))

print("finished high rates")
