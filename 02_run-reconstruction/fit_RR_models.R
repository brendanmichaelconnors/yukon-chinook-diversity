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

# Figure S2 - total run-reconstruction fits
jpeg("./04_figures/figures/figureS2.jpeg", height=4.5, width=7,units="in",res=400,bg = "transparent")

par( mar=c(5,5,1,1) )

#t <- !is.na(rpt$I_t)
t <- rep(TRUE,rpt$nT)
yrF  <- as.factor(rpt$years[t])
I_t <- rpt$I_t[t]*1e-3/exp(rpt$lnqI_s)
E_t <- colSums(exp(rpt$lnRunSize_st))[t]*1e-3
sonarN_t <- colSums(rpt$E_dtg[ ,t,1])*1e-3
ymax <- max(I_t,E_t,sonarN_t,na.rm=TRUE)
counts <- rbind(I_t, (sonarN_t-I_t)); counts[2,25:35] <- sonarN_t[25:35]; counts[is.na(counts)] <-0
b <- barplot(counts, las=1, yaxs="i", xlab="Year",
             ylab="Total border passage (000s)", ylim=c(0,1.1*ymax), col=c(grey(0.45),grey(0.85)) )
#grid()
box(col="grey")

if( is.finite(rpt$sdrpt[1,5]) )
{
  Ese <- filter(rpt$sdrpt,par=="runSize_t")[t, ]
  segments( x0= b, y0=Ese$lCI*1e-3, y1=Ese$uCI*1e-3, col="black", lwd=3 )
}

points(x= b, y= E_t,col="red",cex=0.8)
#points(x= b, y= jtc_data$jtc_bord_pass/1000,col="blue",cex=0.8)

legend( x="topright", bty="n",
        legend=c("Run reconstruction","Mark-recapture","Sonar"),
        pch=c(NA,15,15), pt.cex=c(NA,1.4,1.4), lwd=c(3,NA,NA), col=c("black",grey(0.45),grey(0.85)), lty=c(1,0,0) )
points(x= 29.8, y= 119,col="red",cex=0.8)

axis(1,b[c(1,6,11,16,21,26,31,35)],labels= c("1985","1990","1995","2000", "2005","2010","2015","2019"))
dev.off() 


plotRunSizeMulti()

plotRunSizeMulti <- function( rptFiles=c("./rr_outputs/rpt.base.Rdata",
                                         "./rr_outputs/rpt.oneCor.Rdata",
                                         "./rr_outputs/rpt.fullCor.Rdata") )
  
  
