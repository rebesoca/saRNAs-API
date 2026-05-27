# Script to launch the saRNA Prediction API

library(plumber)

# Get the directory
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)

# Load and run the API
pr <- plumb("plumber.R")

pr$run(host = "0.0.0.0",
       port = 8000,
       swagger = T)