
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))
packages <- c("msm",
              "survival",
              "flexsurv",
              "magrittr",
              "dplyr",
              "data.table",
              "nhm",
              "here",
              "devtools",
              "truncnorm",
              "future",
              "future.apply",
              "Rcpp"
)
installed <- rownames(installed.packages())
to_install <- setdiff(packages, installed)
if (length(to_install) > 0) {
  install.packages(to_install, lib=Sys.getenv("R_LIBS_USER"), repos="https://cloud.r-project.org")
}


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
library(magrittr)

source("fun/models_helpers.R")
source("fun/running_models_fun.R")

rates <- c("low","medium","high")
cores <- 50

# LOW EVENT RATES
system.time(run_models(rates[1],models=c("nhm"),covnames = c("cov1","cov2","cov3"),cores=cores,seeds = 1:50))

print("finished low rates")
# MEDIUM EVENT RATES
gc()
system.time(run_models(rates[2],models=c("nhm"),covnames = c("cov1","cov2","cov3"),cores=cores))
print("finished med rates")

# HIGH EVENT RATES
gc()
system.time(run_models(rates[3],models=c("nhm"),covnames = c("cov1","cov2","cov3"),cores=cores))

print("finished high rates")
