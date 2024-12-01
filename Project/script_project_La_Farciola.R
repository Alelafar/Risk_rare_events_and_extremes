library(evd)
library(ismev)

#saving data and pick Food and Health 

data2022<-read.table(file="Project2021.txt", sep='')
date2022 <- as.Date(rownames(data2022), "%Y%m%d")
data=data2022[,c(2,8)]
rownames(data) <- date2022
head(data)
summary(data)

#plots of data

par(mfrow=c(1,1))
plot(date2022,data$Beer,pch=16,xlab="Day",ylab="Beer",type="l")
plot(date2022,data$Hlth,pch=16,xlab="Day", ylab="Health",type="l")
plot(data$Beer,data$Hlth,pch=16,cex=.45, xlab="Beer", ylab="Health" )

#We want to analyse the extreme falls 
data=-data;
plot(date2022,data$Beer,pch=16,xlab="Day",ylab="Beer",type='l')
plot(date2022,data$Hlth,pch=16,xlab="Day", ylab="Health",type='l')
#day variable
day=c(0:(length(data[,1])-1))

#analysis data$Beer
dat=data$Beer
n=length(dat)

#first approach: fit with stationary models
#Endpoints for diagnostics - need at least 30 points for fit
qu.min <- quantile(dat, 0.5)
qu.max <- quantile(dat,(length(dat)-30)/length(dat))
#Mean residual life plot
mrlplot(dat, tlim=c(qu.min,qu.max))
par(mfrow=c(1,2))
tcplot(dat,tlim=c(qu.min, qu.max))
par(mfrow=c(1,1))
#threshold choice
th=quantile(dat, 0.99)
th
#plot excedeences (useful to see clusters)
exc=dat[which(dat >= th)]
plot(date2022[which(dat >= th)],exc,pch=16, xlab="Day", ylab="Exceedances")
#GPD fitting
fit1<-gpd.fit(dat,threshold=th,npy=253)
gpd.diag(fit1)
par(mfrow=c(1,1))
#CI for parameters
fit1$mle+1.96*fit1$se
fit1$mle-1.96*fit1$se
#anotherway
fit2=fpot(dat,threshold=th,npp=253)
par(mfrow=c(2,2))
plot(fit2)
par(mfrow=c(1,1))
fit2$estimate
#profile
par(mfrow=c(1,2))
plot(profile(fit2))
abline(v=0,col=2,lty=2)
par(mfrow=c(1,1))
#Gumbel model
fit1g<-fpot(dat,threshold=th,shape=0,npp=253)
par(mfrow=c(2,2))
plot(fit1g)
par(mfrow=c(1,1))
fit1g$estimate
#return level
rl1 = th + fit2$estimate[1]/fit2$estimate[2] * ((252.93*fit2$pat)^(fit2$estimate[2]) - 1)
rl10 = th + fit2$estimate[1]/fit2$estimate[2] * ((2529.3*fit2$pat)^(fit2$estimate[2]) - 1)
rl100 = th + fit2$estimate[1]/fit2$estimate[2] * ((25293*fit2$pat)^(fit2$estimate[2]) - 1) 
rl1
rl10
rl100
#Plots
par(mfrow=c(1,3))
gpd.prof(z=fit1, m=1, xlow=3.45, xup=3.9, npy = 253, conf = 0.95)
title("Profile Log-likelihood \n of 1-year Return Level")
abline(v=rl1)  #must coincide with MLE
gpd.prof(z=fit1, m=10, xlow=5.8, xup=8.35, npy = 253, conf = 0.95)
title("Profile Log-likelihood \n of 10-year Return Level")
abline(v=rl10)  #must coincide with MLE
gpd.prof(z=fit1, m=100, xlow=8.3, xup=22, npy = 253, conf = 0.95)
title("Profile Log Likelihood \n of 100-year Return Level")
abline(v=rl100) #must coincide with MLE
par(mfrow=c(1,1))
#validate them
length(which(dat >= rl1))
length(which(dat >= rl10))
length(which(dat >= rl100))


#non-stationary methods

#model with rolling window as threshold (POT approach)

#stability plot of parameters for threshold choice (quantile)
l=21*1;
u=numeric()
exceed=numeric()
alpha0=numeric()
alpha1=numeric()
alpha2=numeric()
xi=numeric()
alpha0_m=numeric()
alpha1_m=numeric()
alpha2_m=numeric()
xi_m=numeric()
alpha0_M=numeric()
alpha1_M=numeric()
alpha2_M=numeric()
xi_M=numeric()
thr=numeric()
for(j in 0:9){
  q=0.9+j*100^(-1)
  thr[j+1]=q
  for (i in (l+1):n){
    u[i-l]=quantile(dat[(i-l):(i-1)],q)
    exceed[i-l]=dat[i]-u[i-l]
  }
  covariates=cbind(dat[(l+1):(n-1)],u[1:(n-l-1)],day[(l+1):(n-1)]) 
  fit3=gpd.fit(dat[(l+2):n],threshold=u[2:(n-l)],ydat=covariates, sigl=c(1,2))
  alpha0[j+1]=fit3$mle[1]
  alpha1[j+1]=fit3$mle[2]
  alpha2[j+1]=fit3$mle[3]
  xi[j+1]=fit3$mle[4]
  alpha0_m[j+1]=fit3$mle[1]-1.96*fit3$se[1]
  alpha0_M[j+1]=fit3$mle[1]+1.96*fit3$se[1]
  alpha1_m[j+1]=fit3$mle[2]-1.96*fit3$se[2]
  alpha1_M[j+1]=fit3$mle[2]+1.96*fit3$se[2]
  alpha2_m[j+1]=fit3$mle[3]-1.96*fit3$se[3]
  alpha2_M[j+1]=fit3$mle[3]+1.96*fit3$se[3]
  xi_m[j+1]=fit3$mle[4]-1.96*fit3$se[4]
  xi_M[j+1]=fit3$mle[4]+1.96*fit3$se[4]
}
par(mfrow=c(2,2))
plot(thr,alpha0,type="b", ylim=c(0,0.4), main="Stability plot for alpha_0",xlab="Quantile")
for(i in 1:10){
  lines(c(thr[i],thr[i]), c(alpha0_m[i],alpha0_M[i]))
}
plot(thr,alpha1,type="b", ylim=c(-0.03,0.1), main="Stability plot for alpha_1",xlab="Quantile")
for(i in 1:10){
  lines(c(thr[i],thr[i]), c(alpha1_m[i],alpha1_M[i]))
}
plot(thr,alpha2,type="b", ylim=c(0.1,0.35), main="Stability plot for alpha_2",xlab="Quantile")
for(i in 1:10){
  lines(c(thr[i],thr[i]), c(alpha2_m[i],alpha2_M[i]))
}
plot(thr,xi,type="b", ylim=c(0,0.25), main="Stability plot for xi",xlab="Quantile")
for(i in 1:10){
  lines(c(thr[i],thr[i]), c(xi_m[i],xi_M[i]))
}
par(mfrow=c(1,1))


#choice of lag (1 month is ~ 21 times)
l=21*3;
u=numeric()
exceed=numeric()
for (i in (l+1):n){
  u[i-l]=quantile(dat[(i-l):(i-1)],0.95)
  exceed[i-l]=dat[i]-u[i-l];
}
plot(date2022[(l+1):n],pmax(exceed[(l+1):n],0),pch=16,type="l")
plot(date2022[(l+1):n],pmax(dat[(l+1):n],0),pch=16,type="l")
#plot threshold
plot(date2022[(l+1):n],u,type='l')
#number of exceedences
n_exc_tot=length(which(dat[(l+1):n] >= u))
n_exc_tot

#fitting with GPD
covariates=cbind(dat[(l+1):(n-1)],u[1:(n-l-1)],day[(l+1):(n-1)]) 
fit3=gpd.fit(dat[(l+2):n],threshold=u[2:(n-l)],ydat=covariates, sigl=c(1,2,3)) 
#trend is very close to 0, so we can ignore day as covariate
fit3=gpd.fit(dat[(l+2):n],threshold=u[2:(n-l)],ydat=covariates, sigl=c(1,2))
gpd.diag(fit3)
# AIC
2*(4+fit3$nllh)

#siglink = exp
fit4=gpd.fit(dat[(l+2):n],threshold=u[2:(n-l)],ydat=covariates, sigl=c(1,2), siglink = exp)
gpd.diag(fit4)
# AIC
2*(4+fit4$nllh)


# annual maxima
maxima=apply(matrix(c(dat,dat[n],dat[n],dat[n],dat[n]),ncol=66), 2, max) # compute maxima (adding final data equals to have an exact ratio)
#plot maxima
plot(c(1950:2015),maxima, pch=16, xlab="Day", ylab = "Annual Maxima")
covariates1=cbind(maxima[1:65])
fit6=gev.fit(maxima[2:66], ydat=covariates1, mul=1, sigl=1, method = "Nelder-Mead")
gev.diag(fit6)
# AIC
2*(5+fit6$nllh)


#monthly maxima
maxima1=apply(matrix(c(dat[32:16663]),ncol=66*12), 2, max) # compute maxima (as before but removing 62 observations, 31 at the beginning e 31 at the end)
#plot maxima
plot(seq(1950,2016,1/12)[-c(793)],maxima1,pch=16, xlab="Day", ylab = "Monthly Maxima")
covariates2=cbind(maxima1[1:791])
fit8=gev.fit(maxima1[2:792], ydat=covariates2, mul=1, sigl = 1,  method = "Nelder-Mead")
gev.diag(fit8)
# AIC
2*(5+fit8$nllh)

#extremal index
par(mfrow=c(2,1)) 
exiplot(dat, tlim=c(0.3,6.4),r=1)
exiplot(exceed[(l+1):n], tlim=c(0,4.4),r=1)
par(mfrow=c(1,1)) 

#return level for PoT model

sigma_hat=cbind(rep(1,dim(covariates)[1]), covariates[,c(1,2)])%*%fit3$mle[1:3]
dat_rescaled= (dat[(l+2):n]-u[2:(n-l)])/sigma_hat
plot(date2022[(l+2):n],dat_rescaled,type="l", xlab="Day", ylab="Scaled Beer")
Dat=dat_rescaled[2:(n-l-1)]
fit3_resc=gpd.fit(Dat, threshold = 0, npy=253)
rl1=  fit3_resc$mle[1]/fit3_resc$mle[2] * ((253*1*fit3_resc$rate)^(fit3_resc$mle[2]) -1)
rl10= fit3_resc$mle[1]/fit3_resc$mle[2] * ((253*10*fit3_resc$rate)^(fit3_resc$mle[2]) -1)
rl100= fit3_resc$mle[1]/fit3_resc$mle[2] * ((253*100*fit3_resc$rate)^(fit3_resc$mle[2]) -1)
rl1
rl10
rl100
par(mfrow=c(1,3)) 
gpd.prof(fit3_resc, xlow=2.98, xup= 3.7, m=1, conf=0.95, npy=253)
abline(v=rl1)  #must coincide with MLE
title("Profile Log-likelihood \n of 1-year Return Level")
gpd.prof(fit3_resc, xlow=5.8, xup= 8.45 , m=10, conf=0.95, npy=253)
abline(v=rl10)  #must coincide with MLE
title("Profile Log-likelihood \n of 10-year Return Level")
gpd.prof(fit3_resc, xlow=8.8, xup= 16.5 , m=100, conf=0.95, npy=253)
abline(v=rl100)  #must coincide with MLE
title("Profile Log-likelihood \n of 100-year Return Level")
par(mfrow=c(1,1)) 
#exponential model
sigma_hat=exp(cbind(rep(1,dim(covariates)[1]), covariates[,c(1,2)])%*%fit4$mle[1:3])
dat_rescaled= (dat[(l+2):n]-u[2:(n-l)])/sigma_hat
Dat=dat_rescaled[2:(n-l-1)]
fit4_resc=gpd.fit(Dat, threshold = 0, npy=253)
rl1=  fit4_resc$mle[1]/fit4_resc$mle[2] * ((253*1*fit4_resc$rate)^(fit4_resc$mle[2]) -1)
rl10= fit4_resc$mle[1]/fit4_resc$mle[2] * ((253*10*fit4_resc$rate)^(fit4_resc$mle[2]) -1)
rl100= fit4_resc$mle[1]/fit4_resc$mle[2] * ((253*100*fit4_resc$rate)^(fit4_resc$mle[2]) -1)
rl1
rl10
rl100
gpd.prof(fit4_resc, xlow=2.98, xup= 3.7, m=1, conf=0.95, npy=253)

#confidence interval
rl_min=c(3.04,6.1,9.5)
rl_max=c(3.64,8.1,14.9)
#transform to real scale for each time t
x1_rs=u[2:(n-l)] + sigma_hat*rl1
x1_min=u[2:(n-l)] + sigma_hat*rl_min[1]
x1_max=u[2:(n-l)] + sigma_hat*rl_max[1]
x10_rs=u[2:(n-l)] + sigma_hat*rl10
x10_min=u[2:(n-l)] + sigma_hat*rl_min[2]
x10_max=u[2:(n-l)] + sigma_hat*rl_max[2]
x100_rs=u[2:(n-l)] + sigma_hat*rl100
x100_min=u[2:(n-l)] + sigma_hat*rl_min[3]
x100_max=u[2:(n-l)] + sigma_hat*rl_max[3]
#validate them
length(which(dat[(l+2):n] >= x1_rs))
length(which(dat[(l+2):n] >= x10_rs))
length(which(dat[(l+2):n] >= x100_rs))
length(which(dat[(l+2):n] >= x1_min))
length(which(dat[(l+2):n] >= x10_min))
length(which(dat[(l+2):n] >= x100_min))
length(which(dat[(l+2):n] >= x1_max))
length(which(dat[(l+2):n] >= x10_max))
length(which(dat[(l+2):n] >= x100_max))

#plot return level
#10-year
plot(date2022[(l+1+7300):(l+1+7600)],pmax(dat[(l+1+7300):(l+1+7600)],0),xlab="Day",ylab="Beer",type="l",col="grey",ylim=c(0,12))
lines(date2022[(l+1+7300):(l+1+7600)], u[7300:7600], col="blue")
lines(date2022[(l+1+7300):(l+1+7600)], x10_rs[7300:7600], col="red")
lines(date2022[(l+1+7300):(l+1+7600)], x10_min[7300:7600] ,lty=3)
lines(date2022[(l+1+7300):(l+1+7600)], x10_max[7300:7600] ,lty=3)
legend("topleft", c("Threshold","10-year return level","Confidence interval"),
       col=c("blue","red","black"), lty=c(1,1,3))

plot(date2022[(l+1+12100):(l+1+12350)],pmax(dat[(l+1+12100):(l+1+12350)],0),xlab="1998",ylab="Beer",type="l",col="grey",ylim=c(0,18))
lines(date2022[(l+1+12100):(l+1+12350)], u[12100:12350], col="blue")
lines(date2022[(l+1+12100):(l+1+12350)], x10_rs[12100:12350], col="red")
lines(date2022[(l+1+12100):(l+1+12350)], x10_min[12100:12350] ,lty=3)
lines(date2022[(l+1+12100):(l+1+12350)], x10_max[12100:12350] ,lty=3)
points(date2022[12268+l+1],dat[12268+l+1],pch=19)
legend("topleft", c("Threshold","10-year return level","Confidence interval"),
       col=c("blue","red","black"), lty=c(1,1,3))


# Return level block maxima
# annual maxima
r10=numeric()
r100=numeric()
l=cbind(rep(1,dim(covariates1)[1]), covariates1[,1])%*%fit6$mle[1:2]
sigma=cbind(rep(1,dim(covariates1)[1]), covariates1[,1])%*%fit6$mle[3:4]
maxima_resc=(maxima[2:66]-l)/sigma
fit6_resc=gev.fit(maxima_resc, method = "Nelder-Mead")
r10_resc=qgev(1-1/10, loc=0 , scale=1 , shape=fit6_resc$mle[3] )
r100_resc=qgev(1-1/100, loc=0 , scale=1 , shape=fit6_resc$mle[3] )
#profile likelihood plots
par(mfrow=c(1,2)) 
gev.prof(fit6_resc, xlow=1.78, xup= 4.45 , m=10, conf=0.95)
abline(v=r10_resc)  #must coincide with MLE
title("Profile Log-likelihood \n of 10-year Return Level")
gev.prof(fit6_resc, xlow=4.1, xup= 14.5 , m=100, conf=0.95)
abline(v=r100_resc)  #must coincide with MLE
title("Profile Log-likelihood \n of 100-year Return Level")
par(mfrow=c(1,1)) 
r_min=c(2,4.6)
r_max=c(3.9,12.9)
r10_rs=l+ sigma*r10_resc
r10_min=l + sigma*r_min[1]
r10_max=l + sigma*r_max[1]
r100_rs=l + sigma*r100_resc
r100_min=l + sigma*r_min[2]
r100_max=l + sigma*r_max[2]
#validate them
length(which(maxima[2:66] >= r10_rs))
length(which(maxima[2:66] >= r100_rs))
length(which(maxima[2:66] >= r10_min))
length(which(maxima[2:66] >= r100_min))
length(which(maxima[2:66] >= r10_max))
length(which(maxima[2:66] >= r100_max))

# monthly maxima
r10=numeric()
r100=numeric()
l=cbind(rep(1,dim(covariates2)[1]), covariates2[,1])%*%fit8$mle[1:2]
sigma=cbind(rep(1,dim(covariates2)[1]), covariates2[,1])%*%fit8$mle[3:4]
maxima1_resc=(maxima1[2:792]-l)/sigma
fit8_resc=gev.fit(maxima1_resc, method = "Nelder-Mead")
r10_resc=qgev(1-1/(10*12), loc=0 , scale=1 , shape=fit8_resc$mle[3] )
r100_resc=qgev(1-1/(100*12), loc=0 , scale=1 , shape=fit8_resc$mle[3] )
#profile likelihood plots
par(mfrow=c(1,2)) 
gev.prof(fit8_resc, xlow=6.35, xup= 9.3 , m=10*12, conf=0.95)
abline(v=r10_resc)  #must coincide with MLE
title("Profile Log-likelihood \n of 10-year Return Level")
gev.prof(fit8_resc, xlow=11, xup= 19.5 , m=100*12, conf=0.95)
abline(v=r100_resc)  #must coincide with MLE
title("Profile Log-likelihood \n of 100-year Return Level")
par(mfrow=c(1,1)) 
r_min=c(6.6,11.7)
r_max=c(9.1,19.1)
r10_rs=l+ sigma*r10_resc
r10_min=l + sigma*r_min[1]
r10_max=l + sigma*r_max[1]
r100_rs=l + sigma*r100_resc
r100_min=l + sigma*r_min[2]
r100_max=l + sigma*r_max[2]
#validate them
length(which(maxima1[2:792] >= r10_rs))
length(which(maxima1[2:792] >= r100_rs))
length(which(maxima1[2:792] >= r10_min))
length(which(maxima1[2:792] >= r100_min))
length(which(maxima1[2:792] >= r10_max))
length(which(maxima1[2:792] >= r100_max))


#save important variable for Beer
datF=dat
fit3F=fit3
exceedF=exceed
uF=u;
covariatesF=covariates


#repeat the previous analysis with dat=data$Health

#analysis data$Health
dat=data$Hlth
n=length(dat)
#choice of lag (1 month is ~ 21 times)
l=21*3;
u=numeric()
exceed=numeric()
for (i in (l+1):n){
  u[i-l]=quantile(dat[(i-l):(i-1)],0.95)
  exceed[i-l]=dat[i]-u[i-l];
}

#fitting with GPD
covariates=cbind(dat[(l+1):(n-1)],u[1:(n-l-1)],day[(l+1):(n-1)])
fit3=gpd.fit(dat[(l+2):n],threshold=u[2:(n-l)],ydat=covariates, sigl=c(1,2))
gpd.diag(fit3)
# AIC
2*(4+fit3$nllh)
fit4=gpd.fit(dat[(l+2):n],threshold=u[2:(n-l)],ydat=covariates, sigl=c(1,2), siglink = exp)
gpd.diag(fit4)
# AIC
2*(4+fit4$nllh)
#extremal index
exiplot(exceed[(l+1):n], tlim=c(0.1,3.5),r=1)
exiplot(dat[(l+1):n], tlim=c(0.3,4.5),r=1)


#save important variable for Health
datH=dat
fit3H=fit4
exceedH=exceed
uH=u;
covariatesH=covariates



#multivariate statistics

#bivariate plot
cor(data$Beer,data$Hlth)

#scaled data 
sigma_hatF=cbind(rep(1,dim(covariatesF)[1]), covariatesF[,c(1,2)])%*%fit3F$mle[1:3]
sigma_hatH=exp(cbind(rep(1,dim(covariates)[1]), covariates[,c(1,2)])%*%fit3H$mle[1:3])
#sigma_hatH=cbind(rep(1,dim(covariatesH)[1]), covariatesH[,c(1,2)])%*%fit3H$mle[1:3]
datF_rescaled= (datF[(l+2):n]-uF[2:(n-l)])/sigma_hatF
datH_rescaled= (datH[(l+2):n]-uH[2:(n-l)])/sigma_hatH
fit0F=gpd.fit(datF_rescaled[2:(n-l-1)], threshold = 0, npy=253)
fit0H=gpd.fit(datH_rescaled[2:(n-l-1)], threshold = 0, npy=253)
plot(datF_rescaled, datH_rescaled, pch=16, cex=.5)
#unit Fréchet transformation
datF_frech=rep(0,length(datF_rescaled))
datF_frech[which(datF_rescaled>0)]=-(log( 1 - fit3F$rate*( 1 + fit3F$mle[4]*
                                   datF_rescaled[which(datF_rescaled>0)])^(-1/fit3F$mle[4]) ) )^(-1)
datF_under=datF_rescaled[-which(datF_rescaled>0)]
prob_underF=numeric()
for(i in 1:length(datF_under)){
  prob_underF[i]=length(datF_under[which(datF_under<= datF_under[i])])/length(datF_rescaled)
}
datF_frech[-which(datF_rescaled>0)]=-(log(prob_underF))^(-1)
datH_under=datH_rescaled[-which(datH_rescaled>0)]
datH_frech=rep(0,length(datF_rescaled))
prob_underH=numeric()
for(i in 1:length(datH_under)){
  prob_underH[i]=length(datH_under[which(datH_under<= datH_under[i])])/length(datH_rescaled)
}
datH_frech[-which(datH_rescaled>0)]=-(log(prob_underH))^(-1)
datH_frech[which(datH_rescaled>0)]=-(log( 1 - fit3H$rate*( 1 + fit3H$mle[4]*
                                   datH_rescaled[which(datH_rescaled>0)])^(-1/fit3H$mle[4]) ) )^(-1)
plot(datF_frech,datH_frech, pch=16, cex=.45, log="xy",xlab="Standardized Beer", ylab="Standardized Health")
xlev=-(log(1-fit3F$rate))^(-1)
ylev=-(log(1-fit3H$rate))^(-1)
abline(v=xlev, lwd=2, col="firebrick1")
abline(h=ylev, lwd=2, col="firebrick1")


#POT approach

bfit1=fbvpot(cbind(datF_rescaled,datH_rescaled), threshold = c(0,0),
             model="log")
bfit1
par(mfrow=c(2,2))
plot(bfit1)
bfit1$param["dep"]+1.96*bfit1$std.err["dep"]
bfit1$param["dep"]-1.96*bfit1$std.err["dep"]

bfit1_bis=fbvpot(cbind(datF_rescaled,datH_rescaled), threshold = c(0,0),
                 model="log", cshape=T)
bfit1_bis
par(mfrow=c(2,2))
plot(bfit1_bis)

#test
bfit1_bis$dev-bfit1$dev  #0.04
qchisq(0.95,1)  #3.84
anova(bfit1,bfit1_bis)

#asymmetric models
bfit2 <- fbvpot(cbind(datF_rescaled,datH_rescaled), threshold = c(0,0), model="alog") 
bfit2
par(mfrow=c(2,2))
plot(bfit2)

#same shape
bbfit2 <- fbvpot(cbind(datF_rescaled,datH_rescaled), threshold = c(0,0), model="alog", cshape=T) 
bbfit2
plot(bbfit2)

#test
bbfit2$dev-bfit2$dev  
qchisq(0.95,1)     
anova(bfit2,bbfit2)

#test against symmetric
bfit1_bis$dev-bbfit2$dev #2.39
qchisq(0.95,2)     #5.99
anova(bbfit2,bfit1_bis)


bfit3 <- fbvpot(cbind(datF_rescaled,datH_rescaled), threshold = c(0,0), model="neglog", cshape=T)
bfit4 <- fbvpot(cbind(datF_rescaled,datH_rescaled), threshold = c(0,0), model="bilog",cshape=T)
bfit5 <- fbvpot(cbind(datF_rescaled,datH_rescaled), threshold = c(0,0), model="ct",cshape=T)
bfit6 <- fbvpot(cbind(datF_rescaled,datH_rescaled), threshold = c(0,0), model="negbilog",cshape=T)
bfit7 <- fbvpot(cbind(datF_rescaled,datH_rescaled), threshold = c(0,0), model="hr",cshape=T)

#bfit8 <- fbvpot(cbind(datF_rescaled,datH_rescaled), threshold = c(0,0), model="aneglog")


bfit1_bis$dev + 2*(length(bfit1_bis$param)-1) #aic1
bbfit2$dev + 2*(length(bbfit2$param)-1) #aic2
bfit3$dev + 2*(length(bfit3$param)-1) #aic3
bfit4$dev + 2*(length(bfit4$param)-1) #aic4
bfit5$dev + 2*(length(bfit5$param)-1) #aic5
bfit6$dev + 2*(length(bfit6$param)-1) #aic6
bfit7$dev + 2*(length(bfit7$param)-1) #aic7
#bfit8$dev + 2*length(bfit8$param) #aic8



#annual maxima
annual.max <- matrix(ncol=2,nrow=66) 
annual.max[,1]=apply(matrix(c(data$Beer,data$Beer[n],data$Beer[n],data$Beer[n],data$Beer[n]),ncol=66), 2, max)
annual.max[,2]=apply(matrix(c(data$Hlth,data$Hlth[n],data$Hlth[n],data$Hlth[n],data$Hlth[n]),ncol=66), 2, max)
par(mfrow=c(2,1)) 
plot(1950:2015,annual.max[,1],ylab="annual maxima Beer",xlab="Years")
plot(1950:2015,annual.max[,2],ylab="monthly maxima Health",xlab="Years")


#One-step estimation - correct standard errors 
fit_biv1 <- fbvevd(annual.max,model="log") 
fit_biv1 
par(mfrow=c(3,2)) 
plot(fit_biv1)
fit_biv1$param["dep"]+1.96*fit_biv1$std.err["dep"]
fit_biv1$param["dep"]-1.96*fit_biv1$std.err["dep"]

#asymmetric logistic
fit_biv2 <- fbvevd(annual.max,model="alog") 
fit_biv2
par(mfrow=c(3,2)) 
plot(fit_biv2)

#test
fit_biv2$dev-fit_biv1$dev  #=0.4
qchisq(0.95,2)             #=5.99
#non possiamo rigettare l'hp di simmetria quindi preferiamo quello che è più semplice come modello

#other models
fit_biv3 <- fbvevd(annual.max,model="neglog")
fit_biv4 <- fbvevd(annual.max,model="bilog")
fit_biv5 <- fbvevd(annual.max,model="ct")
fit_biv6 <- fbvevd(annual.max,model="negbilog")

fit_biv1$dev + 2*length(fit_biv1$param) #aic1
fit_biv2$dev + 2*length(fit_biv2$param) #aic2
fit_biv3$dev + 2*length(fit_biv3$param) #aic3
fit_biv4$dev + 2*length(fit_biv4$param) #aic4
fit_biv5$dev + 2*length(fit_biv5$param) #aic5
fit_biv6$dev + 2*length(fit_biv6$param) #aic6

#common shape (best logistic)
fit_biv1b <- fbvevd(annual.max,model="log", cshape=T) 
fit_biv1b

#test
fit_biv1b$dev-fit_biv1$dev  #=1.93
qchisq(0.95,1)             #=3.84




















