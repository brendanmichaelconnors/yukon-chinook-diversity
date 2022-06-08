#-----------------------------------------------------------------------------#
# calcEnsemble.R        
#
# Calculate ensemble estimates of border passage through time
#-----------------------------------------------------------------------------#

calcEnsemble <- function( mods=c("rpt.base","rpt.oneCor","rpt.fullCor"),
                          relativeWeight=c(1,1,1),
                          nsamp = 1e4 )
{
  I <- length(mods)
  p_i <- standardize(relativeWeight)
  N_i <- ceiling(p_i*nsamp)
  Rse <- list()
  for( i in 1:I )
  {
    load(paste("./02_run-reconstruction/rr_outputs/",mods[i],".Rdata",sep=""))
    Rse[[i]] <- filter(rpt$sdrpt,par=="runSize_st") %>%
                mutate( year = rep(rpt$years,each=rpt$nS,times=1),
                        stock  = rep(rpt$stocks,times=rpt$nT)
                        ) %>%
                select( stock, year, MLE=val, SE=se )
  }
  ens <- data.frame( stock=Rse[[1]]$stock,
                     year=Rse[[1]]$year,
                     "2.5%"=NA, "25%"=NA, "50%"=NA, "75%"=NA, "97.5%"=NA, "SD"=NA )
  colnames(ens)[3:8] <- c("2.5%", "25%", "50%", "75%", "97.5%", "SD")
  for( j in 1:nrow(ens) )
  {
    z <- NULL
    for( i in 1:I )
      z <- c(z,suppressWarnings(rnorm(N_i[i],Rse[[i]]$MLE[j],Rse[[i]]$SE[j])))
    z <- z[!is.na(z)]
    z <- z[z>0]

    ens[j,3:7] <- quantile( z, c(0.025,0.25,0.5,0.75,0.975) )
    ens[j,8] <- (ens[j,7]-ens[j,3])/4
  }

  write.csv(ens,file="./02_run-reconstruction/rr_outputs/ensemble.csv",row.names=FALSE)

}


calcEnsemble()








