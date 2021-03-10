#-----------------------------------------------------------------------------#
# calcEnsemble.R                                                                    #
# Simulation-estimation functions for Yukon River Chinook run reconstruction  #
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








