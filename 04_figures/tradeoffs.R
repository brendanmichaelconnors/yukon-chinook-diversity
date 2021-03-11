# equilibrium tradeoffs -------------------------------------------------------#
alphabeta_df <-  Posteriors.df %>%
  dplyr::select(population,alpha,beta) %>%
  mutate(equilibrium = log(alpha)/beta,
         alpha = alpha) %>%
  group_by(population) %>%
  sample_n(1000,replace=TRUE) %>%
  as.data.frame()

U <- seq(0,1,0.01)

short_posteriors <- plyr::ddply(alphabeta_df,c("population"), function(x) {
  alphas <- x$alpha[1:10000] # 1000 changed from 10000
  betas <- x$beta[1:10000]
  data.frame(alphas,betas)
})

wide.posteriors <- array(NA,dim=c(10000,8,2), dimnames=list(NULL, unique(short_posteriors$population), c("alphas","betas")))

for(i in unique(short_posteriors$population)) {
  xx <- subset(short_posteriors,population==i)
  wide.posteriors[,i,1] <- xx[,2]
  wide.posteriors[,i,2] <- xx[,3]
}

alphas <- 8
betas <- 8

t3 <- array(NA,dim=c(length(U),4,10000))

for(w in 1:10000) {
  draw <- sample(10000,1)
  aa<-seq(1, alphas,1)
  bb <- seq(1, alphas,1)
  for (k in 1: alphas) {
    aa[k] <- log(wide.posteriors[draw,k,1])   # 
    bb[k] <- wide.posteriors[draw,k,2] }
  for (i in 1:length(U)) {
    t1 <- matrix(nrow=length(aa),ncol=4)
    for (j in 1:length(aa)) {
      t1[j,] <- SC.eq(U[i],aa[j],bb[j])
    }
    t3[i,,w] <- apply(t1,2,sum)
    t3[i,3:4,w] <- t3[i,3:4,w]/length(aa)
  }
}

t3.median <- apply(t3,1:2,quantile,probs=c(0.5),na.rm=T)
t3.upper <-apply(t3,1:2,quantile,probs=c(0.9),na.rm=T)
t3.lower <- apply(t3,1:2,quantile,probs=c(0.1),na.rm=T)
t3.median[101,2]<-0
t3.upper[101,2]<-0
t3.lower[101,2]<-0

yyy <- as.numeric(t3.median[,3])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
over.med <- predict(xx)
over.med[over.med <0] =0 
over.med[over.med >1] =1 
over.med[1:8]=0

yyy <- as.numeric(t3.upper[,3])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
over.up <- predict(xx)
over.up[over.up <0] =0 
over.up[over.up >1] =1 
over.up[97:101] =1 
over.up[1:8]=0

yyy <- as.numeric(t3.lower[,3])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
over.low <- predict(xx)
over.low[over.low <0] =0 
over.low[over.low >1] =1 
over.low[97:101] =1 
over.low[1:8]=0

yyy <- as.numeric(t3.median[,4])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
ext.med <- predict(xx)
ext.med[ext.med <0] =0 
ext.med[ext.med >1] =1 
ext.med[1:8]=0

yyy <- as.numeric(t3.upper[,4])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
ext.up <- predict(xx)
ext.up[ext.up <0] =0 
ext.up[ext.up >1] =1 
ext.up[97:101] =1 
ext.up[1:8]=0

yyy <- as.numeric(t3.lower[,4])
xxx <- as.numeric((U*100))
xx <- loess(yyy~xxx)
ext.low <- predict(xx)
ext.low[ext.low <0] =0 
ext.low[ext.low >1] =1 
ext.low[97:101] =1 
ext.low[1:8]=0

t3.upper2 <- as.data.frame(t3.upper)
t3.upper2$strata <- rep(c("upper"),(101))

t3.lower2 <- as.data.frame(t3.lower)
t3.lower2$strata <- rep(c("lower"),(101))

t3.median2 <- as.data.frame(t3.median)
t3.median2$strata <- rep(c("median"),(101))


t3.2 <- cbind(t3.upper2, t3.lower2, t3.median2)

names(t3.2) <- c("upp_V1","upp_V2","upp_V3","upp_V4","up_strata",
                 "low_V1","low_V2","low_V3","low_V4","low_strata",
                 "med_V1","med_V2","med_V3","med_V4","med_strata")

t3.2$U <- U

t3.2$over.med <- over.med
t3.2$over.up <- over.up
t3.2$over.low <- over.low

t3.2$ext.med <- ext.med
t3.2$ext <- ext.up
t3.2$ext.low <- ext.low

t3.3 <- t3.2 %>%
  select(-5,-10,-15)