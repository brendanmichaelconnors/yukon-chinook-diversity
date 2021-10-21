#-----------------------------------------------------------------------------#
# fit_RR_models.R
#
# Fit run-reconstruction models with alternative among population correlation in run-timing
#-----------------------------------------------------------------------------#

source("./02_run-reconstruction/ProcData.R")
source("./02_run-reconstruction/initRR.R")
source("./02_run-reconstruction/plot.R")

# fit models 
rpt <- fitRR(ctlFile = "./02_run-reconstruction/estControlFile.base.txt") 
save(rpt,file="./02_run-reconstruction/rr_outputs/rpt.base.Rdata")
plotAll(rpt = rpt, folder = "./02_run-reconstruction/rr_outputs/base_plots")

rpt <- fitRR(ctlFile = "./02_run-reconstruction/estControlFile.oneCor.txt") 
save(rpt,file="./02_run-reconstruction/rr_outputs/rpt.oneCor.Rdata")
plotAll(rpt = rpt, folder = "./02_run-reconstruction/rr_outputs/oneCor_plots")

rpt <- fitRR(ctlFile = "./02_run-reconstruction/estControlFile.fullCor.txt") 
save(rpt,file="./02_run-reconstruction/rr_outputs/rpt.fullCor.Rdata")
plotAll(rpt = rpt, folder = "./02_run-reconstruction/rr_outputs/fullCor_plots")

# generate ensemble estimates of border passage
source("./02_run-reconstruction/calcEnsemble.R")

plotRunSizeMulti()

plotRunSizeMulti <- function( rptFiles=c("./rr_outputs/rpt.base.Rdata",
                                         "./rr_outputs/rpt.oneCor.Rdata",
                                         "./rr_outputs/rpt.fullCor.Rdata") )
  
  
