#supplementary

#interpolated data length####
mdate=read.csv("https://raw.githubusercontent.com/weecology/PortalData/main/Rodents/moon_dates.csv")

#PB
pbm=mdate%>%filter(!newmoonnumber<278, !newmoonnumber>396)
pbm_na=pbm%>%filter(is.na(censusdate))
length(pbm_na$newmoonnumber)/length(pbm$newmoonnumber)

#PP
ppm=mdate%>%filter(!newmoonnumber<411, !newmoonnumber>526)
ppm_na=ppm%>%filter(is.na(censusdate))
length(ppm_na$newmoonnumber)/length(ppm$newmoonnumber)*100

#intercept-slope covariance plot####
#PP

c1=ggplot(coef_df_PP, aes(x=intercept, y=temp, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("mean temp (lag=1)")+xlab("")+
  ggtitle(expression(italic("C. penicillatus")))+scale_color_manual(values=c(control="#9900cc", removal="#FFcc00"))+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

c2=ggplot(coef_df_PP, aes(x=intercept, y=warm_precip, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("warm precipitation")+
  scale_color_manual(values=c(control="#9900cc", removal="#FFcc00"))+xlab("")+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

c3=ggplot(coef_df_PP, aes(x=intercept, y=cool_precip, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("cool precipitation")+
  scale_color_manual(values=c(control="#9900cc", removal="#FFcc00"))+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

ggarrange(c1,c2,c3, common.legend = T, nrow=3, legend="bottom")

#PB
z1=ggplot(coef_df_PB, aes(x=intercept, y=temp, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("mean temp (lag=1)")+xlab("")+
  ggtitle(expression(italic("C. baileyi")))+scale_color_manual(values=c(control="#9900cc", removal="#FFcc00"))+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

z2=ggplot(coef_df_PB, aes(x=intercept, y=warm_precip, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("warm precipitation")+
  scale_color_manual(values=c(control="#9900cc", removal="#FFcc00"))+xlab("")+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

z3=ggplot(coef_df_PB, aes(x=intercept, y=cool_precip, col=treatment))+geom_point()+theme_classic()+
  geom_smooth(method="lm")+ylab("cool precipitation")+
  scale_color_manual(values=c(control="#9900cc", removal="#FFcc00"))+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

ggarrange(z1,z2,z3, common.legend = T, nrow=3, legend="bottom")

#covariate correlation####

get_temp_warmppt_corr=function(split, model){
  
  analysis_set=analysis(split)
  
  covariates=analysis_set[,3:5]
  
  temp_wprecip=cor.test(analysis_set[,3], analysis_set[,4])$estimate
}

get_temp_coolppt_corr=function(split, model){
  
  analysis_set=analysis(split)
  
  covariates=analysis_set[,3:5]
  
  temp_cprecip=cor.test(analysis_set[,3], analysis_set[,5])$estimate
  
}

get_warmppt_coolppt_corr=function(split, model){
  
  analysis_set=analysis(split)
  
  covariates=analysis_set[,3:5]
  
  wprecip_cprecip=cor.test(analysis_set[,4], analysis_set[,5])$estimate
  
}

#PP
#all data

cor.test(ppcont_covs_dat$meantemp_lag1, ppcont_covs_dat$warm_precip)
cor.test(ppcont_covs_dat$meantemp_lag1, ppcont_covs_dat$cool_precip)
cor.test(ppcont_covs_dat$warm_precip, ppcont_covs_dat$cool_precip)

#at multiple origins
sh1=map(PPcontrol_dat$splits, get_temp_warmppt_corr)
sh2=map(PPcontrol_dat$splits, get_temp_coolppt_corr)
sh3=map(PPcontrol_dat$splits, get_warmppt_coolppt_corr)

sh11=as.vector(unlist(sh1))
sh22=as.vector(unlist(sh2))
sh33=as.vector(unlist(sh3))

par(mfrow=c(1,3))

hist(sh11, main="temperature-warm precipitation", xlab="")
hist(sh22, main ="temperature-cool precipitation", xlab="Pearson's correlation coefficient")
hist(sh33, main="warm-cool precipitation", xlab="")
mtext(expression(italic("C. penicillatus")), side = 3, line = -1.5, outer = TRUE, font=12)

#PB
#all data

cor.test(pbcont_covs_dat$meantemp_lag1, pbcont_covs_dat$warm_precip)
cor.test(pbcont_covs_dat$meantemp_lag1, pbcont_covs_dat$cool_precip)
cor.test(pbcont_covs_dat$warm_precip, pbcont_covs_dat$cool_precip)

zh1=map(PBcontrol_dat$splits, get_temp_warmppt_corr)
zh2=map(PBcontrol_dat$splits, get_temp_coolppt_corr)
zh3=map(PBcontrol_dat$splits, get_warmppt_coolppt_corr)

zh11=as.vector(unlist(zh1))
zh22=as.vector(unlist(zh2))
zh33=as.vector(unlist(zh3))

par(mfrow=c(1,3))

hist(zh11, main="temperature-warm precipitation", xlab="")
hist(zh22, main ="temperature-cool precipitation", xlab="Pearson's correlation coefficient")
hist(zh33, main="warm-cool precipitation", xlab="")
mtext(expression(italic("C. baileyi")), side = 3, line = -1.5, outer = TRUE, font=12)

#parameter covariances####
get_cov_matrix=function(split, model) {
  
  covariance_matrix=invertinfo(model$info.matrix_corrected, silent=T, stopOnError = F)$vcov #compute a covariance matrix 
  #from a given Fisher information matrix by inversion (used info.matrix output of model)
  
}


#PB control model
d1=pmap(list(PBcontrol_dat$splits, PBcontrol_dat$model), get_cov_matrix)

d2=do.call(rbind.data.frame, d1)

parameters=c("Intercept", "beta_1", "beta_12", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,48)

label_rows1=d2%>%cbind(d2, params)%>%select(,7:13)

meantemp_wprecip=label_rows1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecip=label_rows1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecip=label_rows1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

#PB exclosure model
dx1=pmap(list(PBexclosure_dat$splits, PBexclosure_dat$model), get_cov_matrix)

dx2=do.call(rbind.data.frame, dx1)

parameters=c("Intercept", "beta_1", "beta_12", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,48)

label_rows1x=dx2%>%cbind(dx2, params)%>%select(,7:13)

meantemp_wprecipx=label_rows1x%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecipx=label_rows1x%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecipx=label_rows1x%>%select(cool_precip, params)%>%filter(params=="warm_precip")


par(mfrow=c(2,3))
hist(meantemp_wprecip$meantemp_lag1, main ="meantemp-warm precip", xlab="", ylab="control")
hist(meantemp_cprecip$meantemp_lag1, main="meantemp-cool precip", xlab="", ylab="")
hist(cprecip_wprecip$cool_precip, main="cool-warm precip", xlab="", ylab="")
mtext(expression(italic("C. baileyi")), side = 3, line = -1.5, outer = TRUE, font=12)

hist(meantemp_wprecipx$meantemp_lag1, main =NULL, xlab="", ylab="removal")
hist(meantemp_cprecipx$meantemp_lag1, main=NULL, xlab="covariance", ylab="")
hist(cprecip_wprecipx$cool_precip, main=NULL, xlab="", ylab="")
mtext(expression(italic("C. baileyi")), side = 3, line = -1.5, outer = TRUE, font=12)

#PP control model 
f1=pmap(list(PPcontrol_dat$splits, PPcontrol_dat$model), get_cov_matrix)

f2=do.call(rbind.data.frame, f1)

parameters=c("Intercept", "beta_1", "beta_12", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,45)

label_rows2=f2%>%cbind(f2, params)%>%select(,7:13)

PPmeantemp_wprecip=label_rows2%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
PPmeantemp_cprecip=label_rows2%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
PPcprecip_wprecip=label_rows2%>%select(cool_precip, params)%>%filter(params=="warm_precip")

#PP exclosure model
fx1=pmap(list(PPexclosure_dat$splits, PPexclosure_dat$model), get_cov_matrix)

fx2=do.call(rbind.data.frame, fx1)

parameters=c("Intercept", "beta_1", "beta_12", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,45)

label_rows2x=fx2%>%cbind(fx2, params)%>%select(,7:13)

PPmeantemp_wprecipx=label_rows2x%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
PPmeantemp_cprecipx=label_rows2x%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
PPcprecip_wprecipx=label_rows2x%>%select(cool_precip, params)%>%filter(params=="warm_precip")


par(mfrow=c(2,3))
hist(PPmeantemp_wprecip$meantemp_lag1, main ="meantemp-warm precip", xlab="", ylab="control")
hist(PPmeantemp_cprecip$meantemp_lag1, main="meantemp-cool precip", xlab="", ylab="")
hist(PPcprecip_wprecip$cool_precip, main="cool-warm precip", xlab="", ylab="")
mtext(expression(italic("C. penicillatus")), side = 3, line = -1.5, outer = TRUE, font=12)

hist(PPmeantemp_wprecipx$meantemp_lag1, main =NULL, xlab="", ylab="removal")
hist(PPmeantemp_cprecipx$meantemp_lag1, main=NULL, xlab="covariance", ylab="")
hist(PPcprecip_wprecipx$cool_precip, main=NULL, xlab="", ylab="")
mtext(expression(italic("C. baileyi")), side = 3, line = -1.5, outer = TRUE, font=12)

#parameter correlations####

#PP control
pp_c1=pmap(list(PPcontrol_dat$splits, PPcontrol_dat$model), get_cov_matrix)

pp_cd1=map(pp_c1, cov2cor)
pp_cd2=do.call(rbind.data.frame, pp_cd1)

parameters=c("Intercept", "beta_1", "beta_12", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,45)

label_rowsc1=pp_cd2%>%cbind(pp_cd2, params)%>%select(,7:13)

meantemp_wprecipc=label_rowsc1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecipc=label_rowsc1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecipc=label_rowsc1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

#PP exclosure
pp_d1=pmap(list(PPexclosure_dat$splits, PPexclosure_dat$model), get_cov_matrix)

pp_dc1=map(pp_d1, cov2cor)
pp_dc2=do.call(rbind.data.frame, pp_dc1)

parameters=c("Intercept", "beta_1", "beta_12", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,45)

label_rowscx1=pp_dc2%>%cbind(pp_dc2, params)%>%select(,7:13)

meantemp_wprecip2=label_rowscx1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecip2=label_rowscx1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecip2=label_rowscx1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

par(mfrow=c(2,3))
hist(meantemp_wprecipc$meantemp_lag1, main ="meantemp-warm precip", xlab="", ylab="control")
hist(meantemp_cprecipc$meantemp_lag1, main="meantemp-cool precip", xlab="", ylab="")
hist(cprecip_wprecipc$cool_precip, main="cool-warm precip", xlab="", ylab="")
mtext(expression(italic("C. penicillatus")), side = 3, line = -1.5, outer = TRUE, font=12)

hist(meantemp_wprecip2$meantemp_lag1, main =NULL, xlab="", ylab="removal")
hist(meantemp_cprecip2$meantemp_lag1, main=NULL, xlab="correlation", ylab = "")
hist(cprecip_wprecip2$cool_precip, main=NULL, xlab="", ylab="")
mtext(expression(italic("C. penicillatus")), side = 3, line = -1.5, outer = TRUE, font=12)

#PB control
pb_c1=pmap(list(PBcontrol_dat$splits, PBcontrol_dat$model), get_cov_matrix)

pb_cd1=map(pb_c1, cov2cor)
pb_cd2=do.call(rbind.data.frame, pb_cd1)

parameters=c("Intercept", "beta_1", "beta_12", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,48)

label_rowsd1=pb_cd2%>%cbind(pb_cd2, params)%>%select(,7:13)

meantemp_wprecipd=label_rowsd1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecipd=label_rowsd1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecipd=label_rowsd1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

#PB exclosure
pb_d1=pmap(list(PBexclosure_dat$splits, PBexclosure_dat$model), get_cov_matrix)

pb_dc1=map(pb_d1, cov2cor)
pb_dc2=do.call(rbind.data.frame, pb_dc1)

parameters=c("Intercept", "beta_1", "beta_12", "meantemp_lag1", "warm_precip", "cool_precip")
params=rep(parameters,48)

label_rowsdx1=pb_dc2%>%cbind(pb_dc2, params)%>%select(,7:13)

meantemp_wprecipd2=label_rowsdx1%>%select(meantemp_lag1, params)%>%filter(params=="warm_precip")
meantemp_cprecipd2=label_rowsdx1%>%select(meantemp_lag1, params)%>%filter(params=="cool_precip")
cprecip_wprecipd2=label_rowsdx1%>%select(cool_precip, params)%>%filter(params=="warm_precip")

par(mfrow=c(2,3))
hist(meantemp_wprecipd$meantemp_lag1, main ="meantemp-warm precip", xlab="", ylab="control")
hist(meantemp_cprecipd$meantemp_lag1, main="meantemp-cool precip", xlab="", ylab="")
hist(cprecip_wprecipd$cool_precip, main="cool-warm precip", xlab="", ylab="")
mtext(expression(italic("C. baileyi")), side = 3, line = -1.5, outer = TRUE, font=12)

hist(meantemp_wprecipd2$meantemp_lag1, main =NULL, xlab="", ylab="removal")
hist(meantemp_cprecipd2$meantemp_lag1, main=NULL, xlab="correlation", ylab = "")
hist(cprecip_wprecipd2$cool_precip, main=NULL, xlab="", ylab="")
mtext(expression(italic("C. baileyi")), side = 3, line = -1.5, outer = TRUE, font=12)

#models with scaled covariates####
s_covar=covars%>%mutate(meantemps=scale(meantemp_lag1, center = T, scale = T)[, 1],
                        wprecs=scale(warm_precip, center = T, scale = T)[, 1],
                        cprecs=scale(cool_precip, center = T, scale = T)[, 1])

#PP scaled data
ppcont_covs=right_join(s_covar, ppcont_dat)%>%
  select(newmoonnumber, meantemps,
         wprecs, cprecs, PP)%>%rename("abundance"="PP")

ppexcl_covs=right_join(s_covar, ppexcl_dat)%>%
  select(newmoonnumber, meantemps,
         wprecs, cprecs, PP)%>%rename("abundance"="PP")

#select data from Sept 2010- Dec 2019
ppcont_covs_dat=ppcont_covs%>% filter(!newmoonnumber<411, !newmoonnumber>526)
ppexcl_covs_dat=ppexcl_covs%>% filter(!newmoonnumber<411, !newmoonnumber>526)

#interpolate missing data
ppcont_covs_dat$abundance=round_na.interp(ppcont_covs_dat$abundance)
ppexcl_covs_dat$abundance=round_na.interp(ppexcl_covs_dat$abundance)

#PB scaled data
pbcont_covs=right_join(s_covar, pbcont_dat)%>%
  select(newmoonnumber, meantemps,
         wprecs, cprecs, PB)%>%rename("abundance"="PB")

pbexcl_covs=right_join(s_covar, pbexcl_dat)%>%
  select(newmoonnumber, meantemps,
         wprecs, cprecs, PB)%>%rename("abundance"="PB")

#select data from Dec 1999- June 2009
pbcont_covs_dat=pbcont_covs%>%filter(!newmoonnumber<278, !newmoonnumber>396)
pbexcl_covs_dat=pbexcl_covs%>%filter(!newmoonnumber<278, !newmoonnumber>396)

#interpolate missing data
pbcont_covs_dat$abundance=round_na.interp(pbcont_covs_dat$abundance)
pbexcl_covs_dat$abundance=round_na.interp(pbexcl_covs_dat$abundance)

rolling_mod=function(split) {
  
  analysis_set= analysis(split) #get dataframe
  
  fit_model= tsglm(analysis_set[,"abundance"], model = list(past_obs=c(1,12)), 
                   distr = "nbinom", 
                   xreg  = analysis_set[,2:4], 
                   link  = "log")
}

get_preds=function(split, model,y_source) {
  
  analysis_set= analysis(split) #get training data
  assessment_set=assessment(split) #get testing data
  model$ts = y_source$ts #select the initial condition from matching model
  model$response = y_source$response
  preds=predict(model, n.ahead=12, newxreg= assessment_set[,2:4])$pred
}

#model with AR(13)####

rolling_mod=function(split) {
  
  analysis_set= analysis(split) #get dataframe
  
  fit_model= tsglm(analysis_set[,"abundance"], model = list(past_obs=c(1,13)), 
                   distr = "nbinom", 
                   xreg  = analysis_set[,3:5], 
                   link  = "log")
}

#run before building rolling origin bject
n_moons_yr=13
n_yrs=5
n_moons_train=n_moons_yr*n_yrs
n_moons_test=n_moons_yr*1

#run before making forecasts

#PP
m=rep(seq(1:39), each=13) #no.of splits
code1="same"
code2="switched"

#PB
m=rep(seq(1:42), each=13)
code1="same"
code2="switched"
