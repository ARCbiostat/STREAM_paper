
#install.packages(c('arrow', 'fastDummies', 'truncnorm','pROC'))
#install from source file
#install.packages("P:/SWEMED/SWEMED_Research/Caterina/manuscript2/STREAMS_0.0.0.9000.tar.gz", repos = NULL, type = "source")
# pip install -r requirements.txt--> needs to be in the directory
library(reticulate)

py <- import("sys")$executable
system2(py, c("-m", "ensurepip", "--default-pip"))
system2(py, c("-m", "pip", "--version"))
system2(py, c("-m", "pip", "install",
              "numpy>=1.20",
              "torch>=1.12",
              "pyarrow>=7.0",
              "scikit-learn>=1.0",
              "tqdm>=4.0",
              "pandas>=1.3"))

library(STREAMS)



simulation_ready_001 <- readRDS("P:/SWEMED/SWEMED_Research/Caterina/manuscript2/STREAM_paper/Simulation/simulation_results/medium/simulation_ready_001.rds")


cov_vector <- list(
  "0->1" = c("cov1", "cov2", "cov3"),
  "0->2" = c("cov1", "cov2", "cov3"),
  "1->2" = c("cov1", "cov2", "cov3")
)

est <- run_streams(
  data = as.data.frame(simulation_ready_001[[2]]),
  cov_vector = c("cov1", "cov2", "cov3"),
  python = Sys.which("python"),
  pu_args = list(verbose = TRUE)
)


