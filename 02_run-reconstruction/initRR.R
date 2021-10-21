#-----------------------------------------------------------------------------#
# initRR.R
#
# Initialization script for Yukon River Chinook run reconstruction            
#-----------------------------------------------------------------------------#

library(TMB)
library(RColorBrewer)
library(colorRamps)
library(dplyr)
library(reshape2)
library(parallel)
library(here)
library(mvtnorm)
library(ggplot2)
source("02_run-reconstruction/simRR.R")
source("02_run-reconstruction/fitRR.R")
source("02_run-reconstruction/fitSim.R")
source("02_run-reconstruction/tools.R")
compile("02_run-reconstruction/yukonChinookRunRecon.cpp")
