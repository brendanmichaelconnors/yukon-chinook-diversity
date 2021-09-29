#-----------------------------------------------------------------------------#
# procData.R                                                                  #
# Process raw data for Yukon River Chinook run reconstruction                 #
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

rm(list = ls())
library(dplyr)
library(here)
setwd(here())
counts <- read.csv("./02_run-reconstruction/data/borderCounts.csv")
names(counts)[names(counts)=='count_type'] <- 'gear'
counts$gear <- as.character(counts$gear)
gsi <- read.csv('./02_run-reconstruction/data/28Sep2021-update.gsiSamplesAllProbs.csv')
#gsi <- read.csv('./02_run-reconstruction/data/gsiSamplesAllProbs.csv')
#gsi <- read.csv('data/gsiSamples.csv')
stockID <- read.csv('./02_run-reconstruction/data/stockIDs.csv') %>% arrange(plotOrder)
stockID$stockNum <- stockID$plotOrder
gsi <- dplyr::left_join(gsi, stockID, by='region')
gsi$gear <- as.character(gsi$gear)
gsi$gear[gsi$gear=='Fish Wheel'] 		<- 'fishWheel'
gsi$gear[gsi$gear=='Test Fishery'] 		<- 'eagle'
gsi <- subset(gsi, !is.na(julian))
gsi <- subset(gsi, julian <300 & julian >100)
mr <- read.table('./02_run-reconstruction/data/Yukon_chin_border_passage_indices.txt',header=T)
fwDayAdj <- 1
counts$julian[counts$count_type=='fishWheel'] <- counts$julian[counts$count_type=='fishWheel'] +1
gsi$julian_date[gsi$data_label=='YukonRetro'] <- gsi$julian_date[gsi$data_label=='YukonRetro'] +1
stockNames <- stockID$stock
stockRegion <- stockID$region
gsiGear <- c('eagle','fishWheel')
fDay <- 160
lDay <- 285
days <- fDay:lDay
fYear <- min(gsi$year)
lYear <- max(gsi$year)
yrs <- fYear:lYear
mrI_t <- rep(NA,length(yrs))
mrI_t <- mr$mark_recapture
names(mrI_t) <- yrs
n_sdtg <- array(dim=c(length(stockNames),length(fDay:lDay),length(fYear:lYear),length(gsiGear)))
dimnames(n_sdtg) <- list(stockNames=stockNames,julianDay=days,year = yrs,gears=gsiGear)
for(t in 1:length(yrs))
{
	for (g in 1:length(gsiGear))
	{
		tmp <- subset(gsi, year==yrs[t] & 
							gear == gsiGear[g] &
							!is.na(julian) &
							!is.na(region) &
							prob>0) 		
		if(nrow(tmp)==0)
			n_sdtg[,,t,g] <-NA
		smpls <- unique(tmp[,c('julian','sample_num')])
		if(any(table(smpls$sample_num)>1))
		{	
		 errors <- smpls$sample_num[duplicated(smpls$sample_num)]
		  for(err in errors)
		  {
		  	errDays <- table(tmp$julian[tmp$sample_num==err]) %>%
		  			sort(decreasing=TRUE)
		  	tmp$julian[tmp$sample_num==err] <- as.integer(names(errDays)[1])
		  }	
		}
		gsiDat_gt <- tmp
		for(d in 1:length(days))
		{
			tmp <- subset(gsiDat_gt, year== yrs[t] &
								gear == gsiGear[g] &
								julian == days[d] &
								!is.na(julian) &
								!is.na(region) &
								prob>0)
			if(nrow(tmp)==0)
			{
				n_sdtg[,d,t,g] <- NA
			}
			else	
			{
				nSmpls <- length(unique(tmp$sample_num))
				if( sum(tmp$prob) != nSmpls)
					for(smp in unique(tmp$sample_num))
					{
						nProbs <- tmp[tmp$sample_num==smp,]
						if(sum(nProbs$prob) != 1)
						{
							normProbs <- nProbs$prob/sum(nProbs$prob)
							tmp$prob[tmp$sample_num==smp] <- normProbs
						}
					}
				n_s <- dplyr::summarize(group_by(tmp,stockNum),
								 expCounts = sum(prob))
				n_s <- dplyr::left_join(stockID,n_s, by='stockNum')
				n_s$expCounts[is.na(n_s$expCounts)] <-0
				if(round(sum(n_s$expCounts),10) != nSmpls )
					browser(cat('ERROR: sum of normalized GSI probs != sample size'))
				n_sdtg[,d,t,g] <- n_s$expCounts
			}	
		}	
	}
}
idxGear <- c('eagle','fishWheel')
I_dtg <- array(dim=c( length(fDay:lDay),length(fYear:lYear),length(idxGear)))
dimnames(I_dtg) <- list( julianDay=days, year = yrs, gears=idxGear)
for(t in 1:length(yrs))
{
	for (g in 1:length(idxGear))
	{
		tmp <- subset(counts, year==yrs[t] &  gear == idxGear[g] )
		if(nrow(tmp)==0)
			next()
		for(d in 1:length(days))
		{
			tmp <- subset(counts, year== yrs[t] & gear == idxGear[g] & julian == days[d])
			if(nrow(tmp)==0)	
				I_dtg [d,t,g] <- NA
			else
			{
				I_dtg [d,t,g] <- sum(tmp$count, na.rm=T)
			}
		}	
	}
}
chinookYkData <- list(	I_dtg = I_dtg, mrI_t = mrI_t, n_sdtg = n_sdtg, stockNames = stockNames, stockRegion = stockRegion, gears = idxGear, fDay = fDay, lDay = lDay, days = days, fYear = fYear, lYear = lYear, yrs = yrs )
save(chinookYkData, file='./02_run-reconstruction/data/chinookYkData.Rdata')







