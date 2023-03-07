#supplementary

#intercept-slope covariance plot####
#PP

coef_df=as.data.frame(list(ints, b1, b13, temps, warmprec, coolprec))%>%select(treatment, intercept, beta1, beta13, temp,cool_precip, warm_precip)

c1=ggplot(coef_df, aes(x=intercept, y=temp, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("mean temperature(lag=1)")+
  ggtitle("PP")+scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))

c2=ggplot(coef_df, aes(x=intercept, y=warm_precip, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("warm precipitation")+
  scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))


c3=ggplot(coef_df, aes(x=intercept, y=cool_precip, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("cool precipitation")+
  scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))

ggarrange(c1,c2,c3, common.legend = T, nrow=3)

#PB
z1=ggplot(coef_df_PB, aes(x=intercept, y=temp, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("mean temperature(lag=1)")+
  ggtitle("PB")+scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))

z2=ggplot(coef_df_PB, aes(x=intercept, y=warm_precip, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("warm precipitation")+
  scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))

z3=ggplot(coef_df_PB, aes(x=intercept, y=cool_precip, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("cool precipitation")+
  scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))

ggarrange(z1,z2,z3, common.legend = T, nrow=3)

#covariate correlation####
get_temp_warmppt_corr=function(split, model){
  
  analysis_set=analysis(split)
  
  covariates=analysis_set[,5:7]
  
  temp_wprecip=cor.test(analysis_set[,5], analysis_set[,6])$estimate
  #temp_cprecip=cor.test(analysis_set[,5], analysis_set[,7])$estimate
  #wprecip_cprecip=cor.test(analysis_set[,6], analysis_set[,7])$estimate
  
  #correlation=cbind(temp_wprecip, temp_cprecip, wprecip_cprecip)
  }

get_temp_coolppt_corr=function(split, model){
  
  analysis_set=analysis(split)
  
  covariates=analysis_set[,5:7]
  
  temp_cprecip=cor.test(analysis_set[,5], analysis_set[,7])$estimate
  
}

get_warmppt_coolppt_corr=function(split, model){
  
  analysis_set=analysis(split)
  
  covariates=analysis_set[,5:7]
  
  wprecip_cprecip=cor.test(analysis_set[,6], analysis_set[,7])$estimate
  
}

#PP
h1=map(PPcontrol_dat$splits, get_temp_warmppt_corr)
h2=map(PPcontrol_dat$splits, get_temp_coolppt_corr)
h3=map(PPcontrol_dat$splits, get_warmppt_coolppt_corr)

h11=as.vector(unlist(h1))
h22=as.vector(unlist(h2))
h33=as.vector(unlist(h3))

hist(h11, xlab="Pearson's correlation coefficient", main="temperature-warm precipitation")
hist(h22, xlab="Pearson's correlation coefficient", main="temperature-cool precipitation")
hist(h33, xlab="Pearson's correlation coefficient", main="warm-cool precipitation")

par(mfrow=c(1,3))


#PB
s1=map(PBcontrol_dat$splits, get_temp_warmppt_corr)
s2=map(PBcontrol_dat$splits, get_temp_coolppt_corr)
s3=map(PBcontrol_dat$splits, get_warmppt_coolppt_corr)

s11=as.vector(unlist(s1))
s22=as.vector(unlist(s2))
s33=as.vector(unlist(s3))

hist(s11, xlab="Pearson's correlation coefficient", main="temperature-warm precipitation")
hist(s22, xlab="Pearson's correlation coefficient", main="temperature-cool precipitation")
hist(s33, xlab="Pearson's correlation coefficient", main="warm-cool precipitation")

#parameter covariances####
get_cov_matrix=function(split, model) {
  
  covariance_matrix=invertinfo(model$info.matrix, silent=T, stopOnError = F)$vcov #compute a covariance matrix 
  #from a given Fisher information matrix by inversion (used info.matrix output of model)
  
}

#PB control model
d1=pmap(list(PBcontrol_dat$splits, PBcontrol_dat$model), get_cov_matrix)

d2=do.call(rbind.data.frame, d1)

parameters=c("Intercept", "beta_1", "beta_13", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,50)

label_rows1=d2%>%cbind(d2, params)%>%select(,7:13)

meantemp_wprecip=label_rows1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecip=label_rows1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecip=label_rows1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

par(mfrow=c(1,3))
hist(meantemp_wprecip$meantemp_lag1, main ="meantemp-warm precip", xlab="covariance")
hist(meantemp_cprecip$meantemp_lag1, main="meantemp-cool precip", xlab="covariance")
hist(cprecip_wprecip$cool_precip, main="cool-warm precip", xlab="covariance")

#PB exclosure model
d1=pmap(list(PBexclosure_dat$splits, PBexclosure_dat$model), get_cov_matrix)

d2=do.call(rbind.data.frame, d1)

parameters=c("Intercept", "beta_1", "beta_13", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,50)

label_rows1=d2%>%cbind(d2, params)%>%select(,7:13)

meantemp_wprecip=label_rows1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecip=label_rows1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecip=label_rows1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

par(mfrow=c(1,3))
hist(meantemp_wprecip$meantemp_lag1, main ="meantemp-warm precip", xlab="covariance")
hist(meantemp_cprecip$meantemp_lag1, main="meantemp-cool precip", xlab="covariance")
hist(cprecip_wprecip$cool_precip, main="cool-warm precip", xlab="covariance")

#PP control model

d1=pmap(list(PPcontrol_dat$splits, PPcontrol_dat$model), get_cov_matrix)

d2=do.call(rbind.data.frame, d1)

parameters=c("Intercept", "beta_1", "beta_13", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,39)

label_rows1=d2%>%cbind(d2, params)%>%select(,7:13)

meantemp_wprecip=label_rows1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecip=label_rows1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecip=label_rows1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

par(mfrow=c(1,3))
hist(meantemp_wprecip$meantemp_lag1, main ="meantemp-warm precip", xlab="covariance")
hist(meantemp_cprecip$meantemp_lag1, main="meantemp-cool precip", xlab="covariance")
hist(cprecip_wprecip$cool_precip, main="cool-warm precip", xlab="covariance")

#PP exclosure model

d1=pmap(list(PPexclosure_dat$splits, PPexclosure_dat$model), get_cov_matrix)

d2=do.call(rbind.data.frame, d1)

parameters=c("Intercept", "beta_1", "beta_13", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,39)

label_rows1=d2%>%cbind(d2, params)%>%select(,7:13)

meantemp_wprecip=label_rows1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecip=label_rows1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecip=label_rows1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

par(mfrow=c(1,3))
hist(meantemp_wprecip$meantemp_lag1, main ="meantemp-warm precip", xlab="covariance")
hist(meantemp_cprecip$meantemp_lag1, main="meantemp-cool precip", xlab="covariance")
hist(cprecip_wprecip$cool_precip, main="cool-warm precip", xlab="covariance")

#parameter correlations####
#PP control
pp_c1=pmap(list(PPcontrol_dat$splits, PPcontrol_dat$model), get_cov_matrix)

pp_cd1=map(pp_c1, cov2cor)
pp_cd2=do.call(rbind.data.frame, pp_cd1)

parameters=c("Intercept", "beta_1", "beta_13", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,39)

label_rowsc1=pp_cd2%>%cbind(pp_cd2, params)%>%select(,7:13)

meantemp_wprecip=label_rowsc1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecip=label_rowsc1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecip=label_rowsc1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

par(mfrow=c(1,3))
hist(meantemp_wprecip$meantemp_lag1, main ="meantemp-warm precip", xlab="correlation")
hist(meantemp_cprecip$meantemp_lag1, main="meantemp-cool precip", xlab="correlation")
hist(cprecip_wprecip$cool_precip, main="cool-warm precip", xlab="correlation")

#PP exclosure
pp_d1=pmap(list(PPexclosure_dat$splits, PPexclosure_dat$model), get_cov_matrix)

pp_cd1=map(pp_d1, cov2cor)
pp_cd2=do.call(rbind.data.frame, pp_cd1)

parameters=c("Intercept", "beta_1", "beta_13", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,39)

label_rowsc1=pp_cd2%>%cbind(pp_cd2, params)%>%select(,7:13)

meantemp_wprecip=label_rowsc1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecip=label_rowsc1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecip=label_rowsc1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

par(mfrow=c(1,3))
hist(meantemp_wprecip$meantemp_lag1, main ="meantemp-warm precip", xlab="correlation")
hist(meantemp_cprecip$meantemp_lag1, main="meantemp-cool precip", xlab="correlation")
hist(cprecip_wprecip$cool_precip, main="cool-warm precip", xlab="correlation")

#PB control
pb_c1=pmap(list(PBcontrol_dat$splits, PBcontrol_dat$model), get_cov_matrix)

pb_cd1=map(pb_c1, cov2cor)
pb_cd2=do.call(rbind.data.frame, pb_cd1)

parameters=c("Intercept", "beta_1", "beta_13", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,50)

label_rowsc1=pb_cd2%>%cbind(pb_cd2, params)%>%select(,7:13)

meantemp_wprecip=label_rowsc1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecip=label_rowsc1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecip=label_rowsc1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

par(mfrow=c(1,3))
hist(meantemp_wprecip$meantemp_lag1, main ="meantemp-warm precip", xlab="correlation")
hist(meantemp_cprecip$meantemp_lag1, main="meantemp-cool precip", xlab="correlation")
hist(cprecip_wprecip$cool_precip, main="cool-warm precip", xlab="correlation")

#PB exclosure
pb_d1=pmap(list(PBexclosure_dat$splits, PBexclosure_dat$model), get_cov_matrix)

cd1=map(pb_d1, cov2cor)
cd2=do.call(rbind.data.frame, cd1)

parameters=c("Intercept", "beta_1", "beta_13", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,50)

label_rowsc1=cd2%>%cbind(cd2, params)%>%select(,7:13)

meantemp_wprecip=label_rowsc1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecip=label_rowsc1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecip=label_rowsc1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

par(mfrow=c(1,3))
hist(meantemp_wprecip$meantemp_lag1, main ="meantemp-warm precip", xlab="correlation")
hist(meantemp_cprecip$meantemp_lag1, main="meantemp-cool precip", xlab="correlation")
hist(cprecip_wprecip$cool_precip, main="cool-warm precip", xlab="correlation")