#-----------------------------------------------------------------------------#
# fitSim.R  
#
# Simulation-estimation functions for Yukon River Chinook run reconstruction  
#-----------------------------------------------------------------------------#

fitSim <- function( estFolder="mod1",
                    simFolder="sim_N1",
                    nSim=100,
                    nParallelCores=8 )
{
  suppressWarnings(dir.create(paste("simOutputs",simFolder,estFolder,sep="/")))
  mods <- 1:nSim
  cl <- makeCluster(nParallelCores)
  parSapply( cl=cl, X=mods, FUN=fitIndividualSim, simFolder=simFolder, estFolder=estFolder )
  stopCluster(cl)
}

fitIndividualSim <- function( i, simFolder, estFolder )
{
  source("initRR.R")
  ctlFile <- paste(estFolder,"/estControlFile.txt",sep="")
  load(paste("simOutputs/",simFolder,"/sim",i,".Rdata",sep=""))
  rpt <- fitRR( ctlFile=ctlFile, folder=".", simData=sim$simData, saveRun=FALSE )
  save(rpt,file=paste("simOutputs/",simFolder,"/",estFolder,"/fit",i,".Rdata",sep=""))
}

processSimEst <- function( ctrl=NULL,
                           simFolder="sim_N1",
                           estFolder="mod1-uncor",
                           nSim=100 )
{
  if( is.null(ctrl) )
  {
    controlTable  <- .readParFile( "estControlFile.txt" )
    ctrl <- .createList( controlTable )
  }
  converge <- numeric(nSim)
  simDir <- paste('simOutputs/',simFolder,sep="")
  estDir <- paste('simOutputs/',simFolder,"/",estFolder,sep="")
  simFiles <- list.files(simDir, pattern='sim')
  fitFiles <- list.files(estDir, pattern='fit')
  if(length(fitFiles) != length(simFiles))
    cat('WARNING: # sims =',length(simFiles),
                 '| # fits =',length(fitFiles),'\n')
  fitNums <- grep(fitFiles,pattern='fit')
  for( i in fitNums )
  {
    if( is.na(i) ) browser()
    load(paste(simDir,"/sim",i,".Rdata",sep=""))
    load(paste(estDir,"/fit",i,".Rdata",sep=""))
    cvg <- substr(rpt$opt$message,1,4)
    if(length(cvg)>0)
      converge[i] <- cvg %in% c("sing","rela")
  }
  fitNums <- fitNums[converge==1]
  stocks <- ctrl$stocks
  yrs <- ctrl$initYear:ctrl$lastYear
  nReps <- length(fitNums)
  nStocks <- length(stocks)
  nYrs    <- length(yrs)
  om    <- list()
  om$mu_ist <- array(NA, dim=c(nReps,nStocks,nYrs))
  om$N_ist <- array(NA, dim=c(nReps,nStocks,nYrs))
  om$filename <- vector(length=nReps)
  est   <- list()
  est$mu_ist   <- array(NA, dim=c(nReps,nStocks,nYrs))
  est$N_ist    <- array(NA, dim=c(nReps,nStocks,nYrs))
  est$filename <- vector(length=nReps)
  est$resids_its <- array(NA, dim=c(nReps,nYrs,nStocks))
  est$relErr_its <- array(NA, dim=c(nReps,nYrs,nStocks))
  est$muRelErr_its <- array(NA, dim=c(nReps,nYrs,nStocks))
  stats <- list()
  stats$RMSE_is <- array(NA, dim=c(nReps,nStocks))
  stats$MRE_is <- array(NA, dim=c(nReps,nStocks))
  stats$muMRE_is <- array(NA, dim=c(nReps,nStocks))
  stats$CV_is <- array(NA, dim=c(nReps,nStocks))
  stats$stocks <- stocks
  for( i in 1:length(fitNums) )
  {
    load(paste(simDir,"/sim",fitNums[i],".Rdata",sep=""))
    omMu_st <- sim$om$mu_st
    omN_st <- apply(sim$om$N_sdt, c(1,3), sum, na.rm=T)
    om$mu_ist[i,,] <- omMu_st
    om$N_ist[i,,] <- omN_st
    om$filename <- paste("sim",fitNums[i],".Rdata",sep="")
    load(paste(estDir,"/fit",fitNums[i],".Rdata",sep=""))
    estMu_st <- rpt$mu_st
    estN_st <- apply(rpt$N_dst, c(2,3), sum, na.rm=T)
    est$mu_ist[i,,] <- estMu_st
    est$N_ist[i,,] <- estN_st
    est$filename <- paste("fit",fitNums[i],".Rdata",sep="")
    for (s in 1:nStocks)
    {
      muResids <- omMu_st[s,]-estMu_st[s,]
      muRelErr <- muResids/omMu_st[s,]
      resids <- omN_st[s,]-estN_st[s,]
      relErr <- resids/omN_st[s,]
      est$resids_its[i,,s] <- resids
      est$relErr_its[i,,s] <- relErr
      est$muRelErr_its[i,,s] <- muRelErr
      MSE <- mean(resids^2)
      stats$RMSE_is[i,s] <- sqrt(MSE)
      stats$MRE_is[i,s] <- median(relErr)
      stats$muMRE_is[i,s] <- median(muRelErr)
      se <- sd(resids)/sqrt(length(resids))
      stats$CV_is[i,s] <- sd(resids)
    }
  }
  perf <- list(reps = fitNums, om = om, est = est, stats = stats)
  filename <- file.path(estDir,'perf.Rdata')
  save(perf, file=filename)
}



