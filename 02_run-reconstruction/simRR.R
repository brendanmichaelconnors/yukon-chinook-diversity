#-----------------------------------------------------------------------------#
# simRR.R                                                                     
#
# Simulation functions for Yukon River Chinook run reconstruction             
#-----------------------------------------------------------------------------#

simRR <- function( ctlFile = "simCtlFile.txt", folder=NULL,
                   graphics=TRUE)
{ 
  ctl <- .readParFile( ctlFile )
  ctl <- .createList( ctl )
  if ( is.null(folder) )
  {
    stamp <- paste( format(Sys.time(),format="%d%m%Y%H%M%S" ),sep="" )
    folder <- paste ( "simOutputs/sim_",stamp, sep = "" )
  } else folder <- paste ("simOutputs/sim_", folder, sep = "")
  if( !dir.exists(folder) )
      dir.create(folder)
  file.copy(from=ctlFile,to=file.path(folder,"ctlFile.txt"))
  nReps <-  ctl$om$nReps
  for(r in 1:nReps)
  {  
    sim <- chinookOpMod(ctlList = ctl)
    simName <- paste('sim',r,'.Rdata',sep='')
    save(sim,file=file.path(folder,simName))
    if(graphics & r==1)
    {  
      for (yr in c(1985,2008,2016))
        genGraphics(rep=sim, plotYear=yr,
                  graphicsFolder=file.path(folder,'plots'))
    }
  }  

}

chinookOpMod <- function (ctlList = ctl, quiet=TRUE)
{  
  om    <- ctlList$om
  data  <- ctlList$data
  mp    <- ctlList$mp
  fDay <- om$fDay
  lDay <- om$lDay
  day <- seq(from=fDay,to=lDay,by=1)
  nDays <- length(day)
  fYear <- om$fYear
  lYear <- om$lYear
  yrs <- seq(from=fYear,to=lYear,by=1)
  nYrs <- length(yrs)
  stockID <- read.csv(data$stockID)
  nStocks <- nrow(stockID)
  gears <- data$gearNames
  nG    <- length(gears)
  if( !is.null(om$fittedMod) )
    load(paste(om$fittedMod,"/rpt.Rdata",sep=""))
  meanTotN <- om$meanTotN
  if(!is.null(data$propInput))
  {
    load(data$propInput)
    n_sdtg <- chinookYkData$n_sdtg
    n_stg <- apply(n_sdtg,c(1,3,4),sum,na.rm=T)
    n_yrs <- chinookYkData$yrs
    propYrs <- data$propYrs
    propG   <- data$propG
    P_stg   <- array(NA,dim=c(nStocks,
                              length(n_yrs),
                              length(propG)))
    yrIndex <- which(chinookYkData$yrs %in% propYrs)
    for (g in propG)
      for(t in yrIndex)
      {
          n_s <- n_stg[,t,g]
          P_stg[,t,g] <- n_s/sum(n_s)
      } 
    if(propG >1 ) 
      P_st <- apply(P_stg,c(1,2),mean,na.rm=T)
    if(propG ==1)
      P_st <- P_stg[,,propG]
      P_s <- apply(P_st[,yrIndex],1,mean)
     if(sum(P_s) !=1 )
      P_s <- P_s/sum(P_s)
  }  else if(is.null(data$propInput))
    P_s <- 1/nStocks
  N_mu <- log(meanTotN*P_s)
  sigmaN <- om$sigmaN
  lnRunSize_st <- array(NA,dim=c(nStocks,nYrs))
  for(t in 1:nYrs)
    lnRunSize_st[,t] <- N_mu + rnorm(n=nStocks,mean=0,sd=sigmaN)
  if( !is.null(om$fittedMod) )
    lnRunSize_st <- rpt$lnRunSize_st
  runSize_st <- exp(lnRunSize_st)
  mu_st   <- array(NA,c(nStocks,nYrs))
  arrivSD_s <- rep(ctlList$om$arrivSD,nStocks)
  if(!is.null(data$mu_sInput))
  {
    gsi <- read.csv(data$mu_sInput)
    if(data$mu_sInput == "data/gsiSamplesAllProbs.csv")
    {  
      gsi <-  subset(gsi,  sample_type=='gill_eagle'&
                          year>2006) %>% 
              na.omit()
      smry <- dplyr::summarize(group_by(gsi, region,year), n = sum(prob), mean = weighted.mean(x=julian_date, w=prob)
                    )
    }  
    if(data$mu_sInput == "data/gsiSamples.csv")
    {  
      gsi <-  subset(gsi,  gear=='gillnet'&
                           year>2006) %>%
              na.omit()
      smry <- dplyr::summarize(group_by(gsi, region,year),
                    n    = length(unique(sample_num)),
                    mean = mean(julian_date)
                    )
    }  
    meanDay <- dplyr::summarize(group_by(smry,region),
                      avg=mean(mean)
                      )
    meanDay <- dplyr::left_join(meanDay,
                      stockID[,c('region','id')],
                      by='region') %>%
                arrange(id)         
    for(t in 1:nYrs)
      mu_st[,t] <- meanDay$avg
  } else 
  { 
    for(t in 1:nYrs)
      mu_st[,t] <- sample(x=190:210,size=nStocks,replace=T)
  } 
  if( !is.null(om$fittedMod) )
  {
    mu_st <- rpt$mu_st
    if( om$muRW )
    {
      sigma <- matrix( data=0.3, nrow=nStocks, ncol=nStocks )
      diag(sigma) <- 0.035
      for( t in 2:nYrs )
        mu_st[ ,t] <- rmvnorm( n=1, mean=mu_st[ ,t-1], sigma=sigma ) %>%
                      as.numeric()
    }
  }
  N_sdt <- array(NA, dim=c(nStocks,nDays,nYrs))
  rho <- array(NA, dim=c(nStocks,nDays))
  tmpRho <- vector("numeric",length=nDays)
  for(t in 1:nYrs)
  {  
    for( s in 1:nStocks )
    {
        for( d in 1:nDays )
        {
          tmpRho[d] <- exp(-1./(2*arrivSD_s[s])*(day[d]-mu_st[s,t])^2)
        }
        rho[s,] <- tmpRho/sum(tmpRho)
        for( d in 1:nDays )
        {  
            N_sdt[s,d,t] <- runSize_st[s,t]*rho[s,d]
        }    
    }
  }
  N_dt <- apply(X=N_sdt,MARGIN=c(2,3),FUN=sum)
  tauMR <- om$tauMR
  qFW   <- om$q_g[2]
  qMR   <- om$qMR
  sonarIdxYrs <- mp$indexYrs_g1 
  fwIdxYrs <- mp$indexYrs_g2 
  mrIdxYrs <- mp$indexYrs_g3 
  N_dt_obsSonar <- array(NA, c(nDays,nYrs))
  N_dt_obsFW <- array(NA, c(nDays,nYrs))
  N_t_obsMR <- rep(NA,nYrs)
  for (t in 1:nYrs)
  {  
    if(yrs[t] %in% sonarIdxYrs)
    {
      obsErrSonar <- rnorm(n=nDays,mean=0,sd=1)
      N_dt_obsSonar[,t] <- rpois( rep(1,nDays), N_dt[ ,t] )
    }
    if(yrs[t] %in% fwIdxYrs)
    {
      obsErrFw <- rnorm(n=nDays,mean=0,sd=1)
      N_dt_obsFW[,t] <- rpois( rep(1,nDays), qFW*N_dt[ ,t] )
    }  
    if(yrs[t] %in% mrIdxYrs)
    {
      obsErrMR <- rnorm(1,mean=0,sd=1)
      N_t_obsMR[t] <- qMR*sum(N_dt[,t])*exp( tauMR*obsErrMR - tauMR*tauMR/2.)
    }  
  }
  P_sdt <- array(0,c(nStocks,nDays,nYrs))
  for(t in 1:nYrs) 
    for( d in 1:nDays )
      if( N_dt[d,t] )
        P_sdt[,d,t] <- N_sdt[,d,t]/N_dt[d,t]
  E_dtg <- array(NA,dim=c(nDays,nYrs,nG))
  E_dtg[,,1] <- N_dt_obsSonar
  E_dtg[,,2] <- N_dt_obsFW
  I_t <- N_t_obsMR
  n_dtg <- array(NA, dim=c(nDays,nYrs,nG))
  propN_dt <- array(NA, dim=c(nDays,nYrs))
  if(mp$gsiTiming=='perfect')
    for(t in 1:nYrs)
      propN_dt[,t] <- N_dt[,t]/sum(N_dt[,t])
  load("data/chinookYkData.Rdata") 
  nHist_sdtg <- chinookYkData$n_sdtg[ , ,as.character(fYear:lYear), ]
  n_dtg <- apply( nHist_sdtg, 2:4, sum )
  rtList <- calcRunTiming()
  mu_tg <- rtList$mu_tg
  sd_tg <- rtList$sd_tg
  yrs_g <- list( which(yrs %in%sonarIdxYrs), which(yrs %in%fwIdxYrs) )
  for( g in 1:nG )
  {
    for( t in yrs_g[[g]] )
    {
      if( !is.na(mp$gsiN_g[g]) & !is.na(mu_tg[t,g]) )
      {
        p <- dnorm(x=day,mean=mu_tg[t,g],1.5*sd_tg[t,g])
        n_dtg[ ,t,g] <- round( mp$gsiN_g[g] * p )
        x <- mp$gsiN_g[g] - sum(n_dtg[ ,t,g])
        d <- which(day==round(mu_tg[t,g]))
        n_dtg[d,t,g] <- n_dtg[d,t,g] + x
      }
    }
  }
  gsiYrsEagle <- mp$gsiYrs_g1
  gsiYrsFW    <- mp$gsiYrs_g2
  n_sdtg <- array(NA, dim=c(nStocks,nDays,nYrs,nG))
  for(t in 1:nYrs)
    for(d in 1:nDays)
      for(g in 1:nG)
        if(!is.na(n_dtg[d,t,g]))
        {
          if( n_dtg[d,t,g] )
            n_sdtg[,d,t,g]  <- rmultinom(1,n_dtg[d,t,g],P_sdt[ ,d,t])
          else
            n_sdtg[,d,t,g]  <- 0
        }
  simData <- list(n_sdtg = n_sdtg, E_dtg = E_dtg, I_t = I_t, day_d = day )
  om      <- list(N_sdt = N_sdt, mu_st = mu_st, P_sdt = P_sdt, julianDay = day, yrs = yrs, stockNames = stockID$stock, stockRegion = stockID$region)
  report <- list(simData = simData, om = om, ctlList = ctlList)
  return(report)
}

.calcPs <- function( uS,tauP=0,epss )
{
  nStocks <- length( epss )
  xs <- vector( mode="numeric",length=nStocks )
  ps <- vector( mode="numeric",length=nStocks )
  xs <- uS
  xs[ xs==0 ] <- 1.e-6
  xs <- log(xs) + tauP*epss - mean( log(xs)+tauP*epss )
  ps <- exp(xs) / sum( exp(xs) )
}

genGraphics <- function(rep, plotYear,
                        graphicsFolder='plots')
{
  if( !dir.exists(graphicsFolder) )
      dir.create(graphicsFolder)  
  day <- rep$om$julianDay
  stockNames <- as.character(rep$om$stockNames)
  nStocks    <- length(rep$om$stockNames)
  yrs   <- rep$om$yrs
  t     <- which(yrs == plotYear)
  N_sd   <- rep$om$N_sdt[,,t]
  N_d    <- apply(N_sd,2,sum)
  N_d_obs_sonar <- rep$simData$E_dtg[,t,1]/1000
  N_d_obs_fw    <- rep$simData$E_dtg[,t,2]/1000
  n_sdg  <- rep$simData$n_sdtg[,,t,]
  pS_sdg <- array(NA,dim(n_sdg))
  for (g in 1:2)
    for( d in 1:length(day))
      pS_sdg[,d,g] <- n_sdg[,d,g]/sum(n_sdg[,d,g])
  sClrs <- c('#1b9e77','#d95f02','#7570b3','#e7298a',
              '#a6cee3','#e6ab02','#a6761d','#666666')
  plotFile <- paste('N_sdI_dg.pdf',plotYear,'.pdf',sep='')
  pdf(file=paste(graphicsFolder,"/", plotFile, sep="") )
  yMax <- max(max(N_d),max(N_d_obs_sonar),na.rm=T)
  XLIM <- c(170,230)
  plot( x=day, y=N_d,type="n",
        ylim=c(0,yMax),las=1,xlim=XLIM,
        ylab="Number returning (000s)",
        xlab="Julian day")
  rect(xleft=day-.5,xright=day+.5,ybottom=0,ytop=N_d,
       border='grey70',col="gray90",density=100)
  for( s in 1:nStocks )
      lines(x=day,y=N_sd[s,], col=sClrs[s], lwd=1.5)
  points(x=day,y=N_d_obs_sonar, pch=19, cex=0.5)
  points(x=day,y=N_d_obs_fw, pch=19, col="red",cex=0.5)
  legend(x=min(XLIM),y=yMax, bty="n",
          legend=c("SONAR","FISHWHEEL", "Total abundance",
                    stockNames),
          col=c("black","red","gray90",sClrs),
          lty=c(NA,NA,rep(1,length(sClrs)+1)),
          pch=c(19,19,rep(NA,length(sClrs)+1)),
          cex=0.8
          )
  dev.off()
  plotFile <- paste('ObsP_sdg',plotYear,'.pdf',sep='')
  pdf(file=paste(graphicsFolder,"/", plotFile, sep="") )
  yMax <- max(max(N_d_obs_fw),max(N_d_obs_sonar),na.rm=T)
  plot( x=day, y=N_d_obs_sonar,type="n",
        ylim=c(0,yMax),las=1, xlim=XLIM,
        ylab="Obs. number returning (000s)",
        xlab="Julian day")
  rect(xleft=day-.5,xright=day+.5,ybottom=0,ytop=N_d_obs_sonar,
       border='grey70', col="gray90",density=100)
  rect(xleft=day-.5,xright=day+.5,ybottom=0,ytop=N_d_obs_fw,
       border='#fc9272',col="#fcbba1",density=100)
  for( s in 1:nStocks )
      lines(x=day,y=yMax*pS_sdg[s,,2],lty=3, col=sClrs[s])
  for( s in 1:nStocks )
      lines(x=day,y=yMax*pS_sdg[s,,1],lty=1,col=sClrs[s])
  x2Labs <- seq(0,1,length=6)
  axis(side=4,at=seq(from=0, to=yMax,length=length(x2Labs)),
      labels=x2Labs,las=1)
    legend(x=min(XLIM),y=yMax, bty="n",
          legend=c("SONAR","FISHWHEEL",stockNames),
          col=c("grey90","#fcbba1",sClrs),
          lty=rep(1,length(sClrs)+2),
          lwd=c(5,5,rep(1,length(sClrs))),
          cex=0.8
          )
  dev.off()
}  

calcRunTiming <- function( n_dtg=NULL )
{
  if( is.null(n_dtg) )
  {
    load("data/chinookYkData.Rdata")
    n_dtg <- apply( chinookYkData$n_sdtg[ , ,-(1:3), ], 2:4, sum )
  }
  days <- as.numeric(names(n_dtg[ ,1,1]))
  nD <- length(days)
  nT <- dim(n_dtg)[2]
  nG <- dim(n_dtg)[3]
  mu_tg <- matrix( data=NA, nrow=nT, ncol=nG )
  sd_tg <- matrix( data=NA, nrow=nT, ncol=nG )
  for( g in 1:nG )
  {
    for( t in 1:nT )
    {
      z <- NULL
      if( sum(n_dtg[ ,t,g],na.rm=1)>0 )
      {
        for( d in 1:nD ) 
        {
          if( !is.na(n_dtg[d,t,g]) )
            z <- append(z,rep(days[d],n_dtg[d,t,g]))
        }
        mu_tg[t,g] <- mean(z)
        sd_tg[t,g] <- sd(z)
      }
    }
  }
  return( list( mu_tg=mu_tg, sd_tg=sd_tg ) )
}