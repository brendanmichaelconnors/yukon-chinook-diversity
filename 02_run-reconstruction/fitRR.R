#-----------------------------------------------------------------------------#
# fitRR.R                                                                     #
# Estimation functions for Yukon River Chinook run reconstruction             #
#                                                                             #
# Copyright 2019 by Landmark Fisheries Research, Ltd.                         #
#                                                                             #
# This software is provided to Essa Technologies in the hope that it will be  #
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty of      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                        #
#                                                                             #
# ALL INTELLECTUAL PROPERTY REMAINS WITH LANDMARK FISHERIES RESEARCH, LTD.    #
# THIS SOFTWARE MAY NOT BE REDISTRIBUTED, SUBLICENCED, COPIED, OR SHARED      #
# OUTSIDE OF ESSA TECHNOLOGIES WITHOUT THE EXPRESS WRITTEN CONSENT OF         #
# LANDMARK FISHERIES RESEARCH, LTD.                                           #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE    #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#-----------------------------------------------------------------------------#


fitRR <- function( ctlFile="estControlFile.txt",arrivSD=NULL,folder="testMod",simData=NULL,saveRun=TRUE )
{
  source("./02_run-reconstruction/initRR.R")
  controlTable  <- .readParFile( ctlFile )
  ctrl <- .createList( controlTable )
  init  <- ctrl$inits
  load("./02_run-reconstruction/data/chinookYkData.Rdata")
  years  <- ctrl$initYear:ctrl$lastYear
  days   <- chinookYkData$days
  gears  <- chinookYkData$gears
  stocks <- ctrl$stocks
  nT     <- length(years)
  nD     <- length(days)
  nG     <- length(gears)
  nS     <- length(stocks)
  suppressWarnings(dir.create(folder))
  n_sdtg <- chinookYkData$n_sdtg[ , ,as.character(years), ]
  E_dtg  <- chinookYkData$I_dtg[ ,as.character(years), ]
  I_t    <- chinookYkData$mrI_t[as.character(years)]
  for( t in 1:nT )
  {
    for( g in 1:nG )
    {
      E <- E_dtg[ ,t,g]
      if( sum(!is.na(E)) > 0 )  E_dtg[is.na(E),t,g] <- 0
    }
  }
  data <- list( n_sdtg = n_sdtg, E_dtg = E_dtg, I_t = I_t, day_d = days )
  if( !is.null(simData) ) data <- simData
  runSize_st <- matrix( data=init$runSize_s, nrow=nS, ncol=nT )
  initMu_s   <- init$arrivMu_s
  sigma_s    <- init$arrivSD_s
  arrivErr_st <- matrix( data=0, nrow=nS, ncol=nT-1 )
  qE_sg      <- matrix( data=1, nrow=nS, ncol=nG )
  qE_sg[ ,2] <- 0.05
  qI_s       <- rep(0.5,nS)
  errSD_s    <- init$errSD_s
  obsErrSD_g <- rep(0.1,nG)
  cor_ss     <- matrix( data=0, nrow=nS, ncol=nS )
  diag(cor_ss) <- 1
  lnDisp_tg  <- matrix( data=log(1e-4), nrow=nT, ncol=nG )
  lnDisp_tg[ ,2] <- log(0.02)
  lnDisp_tg[17,2] <- log(1.5)
  pars <- list( lnRunSize_st = log(runSize_st),lnArrivMu_s = log(initMu_s),lnArrivSD_s = log(sigma_s),arrivErr_st = arrivErr_st,lnErrSD_s = log(errSD_s),logitCor_ss = logit(cor_ss,lb=-1,ub=1),lnqE_sg = log(qE_sg),lnqI_s = log(qI_s),lnWeightI    = log(150),lnDisp_tg    = lnDisp_tg )
  if(!is.null(arrivSD)) pars$lnArrivSD_s <- rep(log(arrivSD),nS)
  qEmap_sg <- matrix( data=1:(nS*nG), nrow=nS, ncol=nG )
  qEmap_sg[ ,1] <- NA
  qEmap_sg[ ,2] <- ctrl$map$qFishWheel_s
  corMap_ss <- matrix( data=1:(nS*nS), nrow=nS, ncol=nS )
  diag(corMap_ss) <- NA
  for( s in 2:nS )
  {
    if( ctrl$map$corType==0 )
      corMap_ss[s,1:(s-1)] <- NA
    else if( ctrl$map$corType==1 )
      corMap_ss[s,1:(s-1)] <- nS^2+1
  }
  corMap_ss <- mirrorMatrix(corMap_ss)
  dispMap <- NA*lnDisp_tg
  map <- list( lnArrivMu_s = as.factor(ctrl$map$arrivMu_s),lnArrivSD_s = as.factor(ctrl$map$arrivSD_s),lnErrSD_s = as.factor(ctrl$map$errSD_s),logitCor_ss = as.factor(corMap_ss),lnqE_sg = as.factor(qEmap_sg),lnqI_s = as.factor(ctrl$map$qI_s),lnWeightI = as.factor(NA),lnDisp_tg = as.factor(dispMap) )
  compile("./02_run-reconstruction/yukonChinookRunRecon.cpp")
  dyn.load(dynlib("./02_run-reconstruction/yukonChinookRunRecon"))
  obj <- MakeADFun( data = data, parameters = pars,map = map,DLL = "yukonChinookRunRecon",random = NULL )
  low <- obj$par*0-Inf
  upp <- obj$par*0+Inf
  optCtrl <- list(  eval.max = ctrl$maxFunEval,  iter.max = ctrl$maxIterations )
  sink("sink.txt")
  opt <- try( nlminb( start = obj$par,objective = obj$fn,gradient = obj$gr,lower = low,upper = upp,control = optCtrl ) )
  sink()
  rptFE <- obj$report()
  sdrep <- NULL
  errGrad <- NULL
  if( mode(opt)=="character" )
  {
    rpt <- obj$report()
    rpt$opt$convergence <- 1
  }
  else
  {
    par <- data.frame( par  = names(opt$par), val  = opt$par, grad = as.numeric(obj$gr()) )
    sink("sink.txt")
    sdobj <- sdreport( obj )
    sdrpt <- summary( sdobj )
    sink()
    if( mode(sdrpt)!="character" )
    {
      colnames(sdrpt) <- c("val","se")
      sdrpt <- as.data.frame(sdrpt) %>%
             mutate( par = rownames(sdrpt), lCI = val - qnorm(.95)*se, uCI = val + qnorm(.95)*se ) %>%
             dplyr::select( par, val, se, lCI, uCI )
      errGrad <- filter( par, par=="arrivErr_st" )
      errGrad_st <- matrix( errGrad$grad, nrow=nS, ncol=nT-1 )
    }
    rpt <- obj$report()
    rpt$opt    <- opt
    rpt$years  <- years
    rpt$gears  <- gears
    rpt$stocks <- stocks
    rpt$par    <- par
    rpt$sdrpt  <- sdrpt
    rpt$errGrad_st <- errGrad_st
    nObs <- sum(!is.na(data$n_sdtg[1, , , ])) + sum(!is.na(data$E_dtg)) +  length(data$I_t)
    nPar <- length(opt$par)
    nll  <- opt$objective
    aic  <- 2*nPar + 2*nll + 2*nPar*(nPar+1)/(nObs-nPar-1)
    rpt$aic <- aic
    if( saveRun )
    {
      save( rpt, file=paste(folder,"/rpt.Rdata",sep="") )
      system( paste("cp ",ctlFile," ",folder,"/estControlFile.txt",sep="") )
    }

  }

  rpt

}