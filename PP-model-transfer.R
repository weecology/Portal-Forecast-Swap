#Script to run model transfer analyses: PP case study###

#PACKAGES####
require(tscount)
require(portalr)
require(dplyr)
require(vctrs)
require(rsample)
require(lubridate)
require(portalcasting)
require(ggplot2)
require(tidymodels)
require(purrr)
require(yardstick)
require(Metrics)
require(tidyr)
require(ggpubr)
require(tidyverse)
require(overlapping)

source("model-transfer-functions.R")

#DATA MANIPULATION####

#rodent data
use_default_data_path("D:\\Dropbox (UFL)\\PhD-stuff\\Portal-Forecast-Swap")

rodent_data=summarize_rodent_data(
  path = get_default_data_path(),
  clean = TRUE,
  level="Treatment",
  type = "Rodents",
  plots = "Longterm",
  unknowns = FALSE,
  shape = "crosstab",
  time = "newmoon",
  output = "abundance",
  na_drop = FALSE,
  zero_drop = FALSE,
  min_traps = 1,
  min_plots = 1,
  effort = TRUE,
  download_if_missing = TRUE,
  quiet = FALSE
)

ppcont_dat=rodent_data%>%
  select(newmoonnumber,treatment, PP)%>%
  replace_na(list(treatment='control'))%>%
  filter(treatment=="control")

ppexcl_dat=rodent_data%>%
  select(newmoonnumber,treatment, PP)%>%
  replace_na(list(treatment='exclosure'))%>%
  filter(treatment=="exclosure")

#covariates data

covars=weather(level="newmoon", fill=TRUE, horizon=365, path=get_default_data_path())%>%
  select(newmoonnumber,meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)%>%
  mutate(meantemp_lag1=lag(meantemp,order_by=newmoonnumber))

ppcont_covs=right_join(covars,ppcont_dat)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PP)%>%rename("abundance"="PP")

ppexcl_covs=right_join(covars,ppexcl_dat)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PP)%>%rename("abundance"="PP")

#select data from Sept 2010- Dec 2019
ppcont_covs_dat=ppcont_covs%>% filter(!newmoonnumber<411, !newmoonnumber>526)
ppexcl_covs_dat=ppexcl_covs%>% filter(!newmoonnumber<411, !newmoonnumber>526)

#interpolate missing data
ppcont_covs_dat$abundance=round_na.interp(ppcont_covs_dat$abundance)
ppexcl_covs_dat$abundance=round_na.interp(ppexcl_covs_dat$abundance)

#rolling origin object for analysis####
n_moons_yr=12
n_yrs=5
n_moons_train=n_moons_yr*n_yrs
n_moons_test=n_moons_yr*1

PPcontrol_dat <- 
  rolling_origin(
    data       = ppcont_covs_dat, #all PP control data (2010-2019)
    initial    = n_moons_train, #samples used for modelling (training)
    assess     = n_moons_test, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed; 
  ) 

PPexclosure_dat <- 
  rolling_origin(
    data       = ppexcl_covs_dat, #all PP exclosure data (2010-2019)
    initial    = n_moons_train, #samples used for modelling (training)
    assess     = n_moons_test, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

#DATA ANALYSES####

#Model fitting####
#fitting matched and mismatched models

#add column for model(same data and model)
PPcontrol_dat$model=map(PPcontrol_dat$splits, rolling_mod)
PPexclosure_dat$model=map(PPexclosure_dat$splits, rolling_mod)

#matched and mismatched model predictions

#add column for model predictions (same data and model)
PPcontrol_dat$preds_same=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$model, PPcontrol_dat$model), get_preds)
PPexclosure_dat$preds_same=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$model, PPexclosure_dat$model), get_preds)

#add column for model predictions from switched model
PPcontrol_dat$preds_switch=pmap(list(PPcontrol_dat$splits, PPexclosure_dat$model, PPcontrol_dat$model), get_preds)
PPexclosure_dat$preds_switch=pmap(list(PPexclosure_dat$splits, PPcontrol_dat$model, PPexclosure_dat$model), get_preds)

#matched and mismatched model evaluations

#add column for model evals from same model
PPcontrol_dat$evals_same=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$preds_same), mod_evals_same)
PPexclosure_dat$evals_same=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$preds_same), mod_evals_same)

#add column for model evals from switched model
PPcontrol_dat$evals_switch=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$preds_switch),mod_evals_switch)
PPexclosure_dat$evals_switch=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$preds_switch),mod_evals_switch)

#Parameter comparison####

ppcontrol_int_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(1)
ppexclosure_int_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(1)

ppcontrol_b1_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(2)
ppexclosure_b1_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(2)

ppcontrol_b12_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(3)
ppexclosure_b12_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(3)

ppcontrol_temp1_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(4)
ppexclosure_temp1_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(4)

ppcontrol_warmprec_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(5)
ppexclosure_warmprec_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(5)

ppcontrol_coolprec_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(6)
ppexclosure_coolprec_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(6)

warmprecpp=cbind(ppcontrol_warmprec_coefs, ppexclosure_warmprec_coefs)%>%as.data.frame%>%
  rename("control"="ppcontrol_warmprec_coefs", "removal"="ppexclosure_warmprec_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "warm_precip")

coolprecpp=cbind(ppcontrol_coolprec_coefs, ppexclosure_coolprec_coefs)%>%as.data.frame%>%
  rename("control"="ppcontrol_coolprec_coefs", "removal"="ppexclosure_coolprec_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "cool_precip")

tempspp=cbind(ppcontrol_temp1_coefs, ppexclosure_temp1_coefs)%>%as.data.frame%>%
  rename("control"="ppcontrol_temp1_coefs", "removal"="ppexclosure_temp1_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "temp")

intspp=cbind(ppcontrol_int_coefs, ppexclosure_int_coefs)%>%as.data.frame%>%
  rename("control"="ppcontrol_int_coefs", "removal"="ppexclosure_int_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "intercept")

b1pp=cbind(ppcontrol_b1_coefs, ppexclosure_b1_coefs)%>%as.data.frame%>%
  rename("control"="ppcontrol_b1_coefs", "removal"="ppexclosure_b1_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "beta1")

b12pp=cbind(ppcontrol_b12_coefs, ppexclosure_b12_coefs)%>%as.data.frame%>%
  rename("control"="ppcontrol_b12_coefs", "removal"="ppexclosure_b12_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "beta12")

coef_df_PP=as.data.frame(list(intspp, b1pp, b12pp, tempspp, warmprecpp, coolprecpp))%>%
  select(treatment, intercept, beta1, beta12, temp,cool_precip, warm_precip)

ppint_cont=coef_df_PP%>%filter(treatment=="control")%>%select(intercept)%>%rename("control"="intercept")
ppint_excl=coef_df_PP%>%filter(treatment=="removal")%>%select(intercept)%>%rename("removal"="intercept")

ppar1_cont=coef_df_PP%>%filter(treatment=="control")%>%select(beta1)%>%rename("control"="beta1")
ppar1_excl=coef_df_PP%>%filter(treatment=="removal")%>%select(beta1)%>%rename("removal"="beta1")

ppar12_cont=coef_df_PP%>%filter(treatment=="control")%>%select(beta12)%>%rename("control"="beta12")
ppar12_excl=coef_df_PP%>%filter(treatment=="removal")%>%select(beta12)%>%rename("removal"="beta12")

pptemp_cont=coef_df_PP%>%filter(treatment=="control")%>%select(temp)%>%rename("control"="temp")
pptemp_excl=coef_df_PP%>%filter(treatment=="removal")%>%select(temp)%>%rename("removal"="temp")

ppcprec_cont=coef_df_PP%>%filter(treatment=="control")%>%select(cool_precip)%>%rename("control"="cool_precip")
ppcprec_excl=coef_df_PP%>%filter(treatment=="removal")%>%select(cool_precip)%>%rename("removal"="cool_precip")

ppwprec_cont=coef_df_PP%>%filter(treatment=="control")%>%select(warm_precip)%>%rename("control"="warm_precip")
ppwprec_excl=coef_df_PP%>%filter(treatment=="removal")%>%select(warm_precip)%>%rename("removal"="warm_precip")

ppo1=cbind(ppint_cont, ppint_excl)
ppo2=cbind(ppar1_cont, ppar1_excl)
ppo3=cbind(ppar12_cont, ppar12_excl)
ppo4=cbind(pptemp_cont, pptemp_excl)
ppo5=cbind(ppcprec_cont, ppcprec_excl)
ppo6=cbind(ppwprec_cont, ppwprec_excl)

##degree of overlap####
overlap(ppo1, plot=T)
overlap(ppo2, plot=T)
overlap(ppo3, plot=T)
overlap(ppo4, plot=T)
overlap(ppo5, plot=T)
overlap(ppo6, plot=T)

##directional shift####

#get est and se for z score calculation

PPcontrol_dat$coef=map(PPcontrol_dat$model, get_coef)
PPexclosure_dat$coef=map(PPexclosure_dat$model, get_coef)

#select estimates for each parameter

ppcont_int_est=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(1)
ppcont_b1_est=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(2)
ppcont_b12_est=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(3)
ppcont_temp_est=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(4)
ppcont_wprec_est=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(5)
ppcont_cprec_est=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(6)

ppexcl_int_est=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(1)
ppexcl_b1_est=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(2)
ppexcl_b12_est=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(3)
ppexcl_temp_est=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(4)
ppexcl_wprec_est=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(5)
ppexcl_cprec_est=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(6)

#select standard errors for each parameter
ppcont_int_se=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(8)
ppcont_b1_se=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(9)
ppcont_b12_se=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(10)
ppcont_temp_se=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(11)
ppcont_wprec_se=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(12)
ppcont_cprec_se=PPcontrol_dat$coef%>%map(unlist)%>%map_dbl(13)

ppexcl_int_se=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(8)
ppexcl_b1_se=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(9)
ppexcl_b12_se=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(10)
ppexcl_temp_se=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(11)
ppexcl_wprec_se=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(12)
ppexcl_cprec_se=PPexclosure_dat$coef%>%map(unlist)%>%map_dbl(13)

#combine dataframes

ppint_coef=cbind(ppcont_int_est, ppexcl_int_est, ppcont_int_se, ppexcl_int_se)%>%
  as.data.frame%>%rename("B_c"="ppcont_int_est", "B_r"="ppexcl_int_est",
                         "SE_c"="ppcont_int_se", "SE_r"="ppexcl_int_se")%>%
  mutate(parameter="intercept", 
         z_score=(B_c - B_r) / sqrt((SE_c + SE_r)^2),
         pvalue=2 * pnorm(z_score, lower.tail = FALSE),
         adj_p=p.adjust(pvalue, method="fdr"),
         B_shift= B_c - B_r,
         significance=case_when(adj_p<0.05 & B_shift<0 ~ "SN", #SN==significant negative shift
                                adj_p<0.05 & B_shift>0 ~ "SP", #SP==significant positive shift
                                adj_p>0.05~ "NS")) #NS==not significant

ppb1_coef=cbind(ppcont_b1_est, ppexcl_b1_est, ppcont_b1_se, ppexcl_b1_se)%>%
  as.data.frame%>%rename("B_c"="ppcont_b1_est", "B_r"="ppexcl_b1_est",
                         "SE_c"="ppcont_b1_se", "SE_r"="ppexcl_b1_se")%>%
  mutate(parameter="beta1", 
         z_score=(B_c - B_r) / sqrt((SE_c + SE_r)^2),
         pvalue=2 * pnorm(z_score, lower.tail = FALSE),
         adj_p=p.adjust(pvalue, method="fdr"),
         B_shift= B_c - B_r,
         significance=case_when(adj_p<0.05 & B_shift<0 ~ "SN", #SN==significant negative shift
                                adj_p<0.05 & B_shift>0 ~ "SP", #SP==significant positive shift
                                adj_p>0.05~ "NS")) #NS==not significant

ppb12_coef=cbind(ppcont_b12_est, ppexcl_b12_est, ppcont_b12_se, ppexcl_b12_se)%>%
  as.data.frame%>%rename("B_c"="ppcont_b12_est", "B_r"="ppexcl_b12_est",
                         "SE_c"="ppcont_b12_se", "SE_r"="ppexcl_b12_se")%>%
  mutate(parameter="beta12", 
         z_score=(B_c - B_r) / sqrt((SE_c + SE_r)^2),
         pvalue=2 * pnorm(z_score, lower.tail = FALSE),
         adj_p=p.adjust(pvalue, method="fdr"),
         B_shift= B_c - B_r,
         significance=case_when(adj_p<0.05 & B_shift<0 ~ "SN", #SN==significant negative shift
                                adj_p<0.05 & B_shift>0 ~ "SP", #SP==significant positive shift
                                adj_p>0.05~ "NS")) #NS==not significant

pptemp_coef=cbind(ppcont_temp_est, ppexcl_temp_est, ppcont_temp_se, ppexcl_temp_se)%>%
  as.data.frame%>%rename("B_c"="ppcont_temp_est", "B_r"="ppexcl_temp_est",
                         "SE_c"="ppcont_temp_se", "SE_r"="ppexcl_temp_se")%>%
  mutate(parameter="temp", 
         z_score=(B_c - B_r) / sqrt((SE_c + SE_r)^2),
         pvalue=2 * pnorm(z_score, lower.tail = FALSE),
         adj_p=p.adjust(pvalue, method="fdr"),
         B_shift= B_c - B_r,
         significance=case_when(adj_p<0.05 & B_shift<0 ~ "SN", #SN==significant negative shift
                                adj_p<0.05 & B_shift>0 ~ "SP", #SP==significant positive shift
                                adj_p>0.05~ "NS")) #NS==not significant

ppwprec_coef=cbind(ppcont_wprec_est, ppexcl_wprec_est, ppcont_wprec_se, ppexcl_wprec_se)%>%
  as.data.frame%>%rename("B_c"="ppcont_wprec_est", "B_r"="ppexcl_wprec_est",
                         "SE_c"="ppcont_wprec_se", "SE_r"="ppexcl_wprec_se")%>%
  mutate(parameter="warm precip", 
         z_score=(B_c - B_r) / sqrt((SE_c + SE_r)^2),
         pvalue=2 * pnorm(z_score, lower.tail = FALSE),
         adj_p=p.adjust(pvalue, method="fdr"),
         B_shift= B_c - B_r,
         significance=case_when(adj_p<0.05 & B_shift<0 ~ "SN", #SN==significant negative shift
                                adj_p<0.05 & B_shift>0 ~ "SP", #SP==significant positive shift
                                adj_p>0.05~ "NS")) #NS==not significant

ppcprec_coef=cbind(ppcont_cprec_est, ppexcl_cprec_est, ppcont_cprec_se, ppexcl_cprec_se)%>%
  as.data.frame%>%rename("B_c"="ppcont_cprec_est", "B_r"="ppexcl_cprec_est",
                         "SE_c"="ppcont_cprec_se", "SE_r"="ppexcl_cprec_se")%>%
  mutate(parameter="cool precip", 
         z_score=(B_c - B_r) / sqrt((SE_c + SE_r)^2),
         pvalue=2 * pnorm(z_score, lower.tail = FALSE),
         adj_p=p.adjust(pvalue, method="fdr"),
         B_shift= B_c - B_r,
         significance=case_when(adj_p<0.05 & B_shift<0 ~ "SN", #SN==significant negative shift
                                adj_p<0.05 & B_shift>0 ~ "SP", #SP==significant positive shift
                                adj_p>0.05~ "NS")) #NS==not significant

#combine dataframes
pp_sig=rbind(ppint_coef, ppb1_coef, ppb12_coef, pptemp_coef, ppwprec_coef, ppcprec_coef)

#Evaluate model transferability####

###forecasts####

#control-control
ppcont_preds_same=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$preds_same), get_dat_same)

#control dat-exclosure mod
ppcont_preds_switch=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$preds_switch), get_dat_switch)

#exclosure
ppexcl_preds_same=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$preds_same), get_dat_same)

#exclosure dat-controlmod
ppexcl_preds_switch=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$preds_switch), get_dat_switch)

m=rep(seq(1:45), each=12) #no.of splits
code1="same"
code2="switched"

PPpreds_cont_same=do.call(rbind.data.frame, ppcont_preds_same)
PPpreds_cont_same=cbind(PPpreds_cont_same, m, code1)
PPpreds_cont_switch=do.call(rbind.data.frame, ppcont_preds_switch)
PPpreds_cont_switch=cbind(PPpreds_cont_switch, m, code2)

PPpreds_excl_same=do.call(rbind.data.frame, ppexcl_preds_same)
PPpreds_excl_same=cbind(PPpreds_excl_same, m, code1)
PPpreds_excl_switch=do.call(rbind.data.frame, ppexcl_preds_switch)
PPpreds_excl_switch=cbind(PPpreds_excl_switch, m, code2)

###RMSE####

ppcont_evals1_diff=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_same, PPcontrol_dat$evals_switch ,PPcontrol_dat$id), get_evals1_diff)
ppevals1_cont_diff=do.call(rbind.data.frame, ppcont_evals1_diff)%>%mutate(plot="control")

ppexcl_evals1_diff=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_same, PPexclosure_dat$evals_switch ,PPexclosure_dat$id), get_evals1_diff)
ppevals1_excl_diff=do.call(rbind.data.frame, ppexcl_evals1_diff)%>%mutate(plot="removal")

#h=1
ppevals1=rbind(ppevals1_cont_diff, ppevals1_excl_diff)
ppevals1$newmoon=as.integer(ppevals1$newmoon)
ppevals1$h=as.integer(ppevals1$h)
ppevals1$score_same=as.numeric(ppevals1$score_same)
ppevals1$score_switch=as.numeric(ppevals1$score_switch)
ppevals1$score_diff=as.numeric(ppevals1$score_diff)

#h=6
ppcont_evals6_diff=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_same, PPcontrol_dat$evals_switch ,PPcontrol_dat$id), get_evals6_diff)
ppevals6_cont_diff=do.call(rbind.data.frame, ppcont_evals6_diff)%>%mutate(plot="control")

ppexcl_evals6_diff=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_same, PPexclosure_dat$evals_switch ,PPexclosure_dat$id), get_evals6_diff)
ppevals6_excl_diff=do.call(rbind.data.frame, ppexcl_evals6_diff)%>%mutate(plot="removal")

ppevals6=rbind(ppevals6_cont_diff, ppevals6_excl_diff)
ppevals6$newmoon=as.integer(ppevals6$newmoon)
ppevals6$h=as.integer(ppevals6$h)
ppevals6$score_same=as.numeric(ppevals6$score_same)
ppevals6$score_switch=as.numeric(ppevals6$score_switch)
ppevals6$score_diff=as.numeric(ppevals6$score_diff)

#h=12
ppcont_evals12_diff=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_same, PPcontrol_dat$evals_switch ,PPcontrol_dat$id), get_evals12_diff)
ppevals12_cont_diff=do.call(rbind.data.frame, ppcont_evals12_diff)%>%mutate(plot="control")

ppexcl_evals12_diff=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_same, PPexclosure_dat$evals_switch ,PPexclosure_dat$id), get_evals12_diff)
ppevals12_excl_diff=do.call(rbind.data.frame, ppexcl_evals12_diff)%>%mutate(plot="removal")

ppevals12=rbind(ppevals12_cont_diff, ppevals12_excl_diff)
ppevals12$newmoon=as.integer(ppevals12$newmoon)
ppevals12$h=as.integer(ppevals12$h)
ppevals12$score_same=as.numeric(ppevals12$score_same)
ppevals12$score_switch=as.numeric(ppevals12$score_switch)
ppevals12$score_diff=as.numeric(ppevals12$score_diff)

###Brier scores####

#control

#combine predictions on same and switched models for control data
pp_preds_control=left_join(PPpreds_cont_same, PPpreds_cont_switch, by=c("moon", "holdout", "m"))

#calculate brier score for same models
pp_brier_cont1=scoring(pred=pp_preds_control$preds_same, response=pp_preds_control$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine same model predictions and brier score dataframes
pp_brier_cont_same=cbind(pp_preds_control, pp_brier_cont1)%>%rename(quadratic_same=quadratic)

#calculate brier score for switched models
pp_brier_cont2=scoring(pred=pp_preds_control$preds_switch, response=pp_preds_control$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine switched model predictions and brier score dataframes
pp_brier_cont_switch=cbind(pp_preds_control, pp_brier_cont2)%>%rename(quadratic_switch=quadratic)

#calculate brier score differences
pp_brier_control=left_join(pp_brier_cont_same, pp_brier_cont_switch)%>%mutate(treatment="control",brier_diff=quadratic_same-quadratic_switch)

#select relevant rows with horizons 1,6 12
pp_brier_control1=pp_brier_control%>%group_by(m)%>%slice(which(row_number()%in%c(1,6,12)))

#subdivide into each horizon for easier plotting
ppb1=pp_brier_control1%>%group_by(m)%>%slice_head(n=1)%>%mutate(horizon="1")
ppb6=pp_brier_control1%>%group_by(m)%>%slice(which(row_number()==2))%>%mutate(horizon="6")
ppb12=pp_brier_control1%>%group_by(m)%>%slice(which(row_number()==3))%>%mutate(horizon="12")

#combine all 3 horizons
ppb=rbind(ppb1,ppb6,ppb12)

#exclosure###

#combine predictions on same and switched models for exclosure data
pp_preds_exclosure=left_join(PPpreds_excl_same, PPpreds_excl_switch, by=c("moon", "holdout", "m"))

#calculate brier score for same models
pp_brier_excl1=scoring(pred=pp_preds_exclosure$preds_same, response=pp_preds_exclosure$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine same model predictions and brier score dataframes
pp_brier_excl_same=cbind(pp_preds_exclosure, pp_brier_excl1)%>%rename(quadratic_same=quadratic)

#calculate brier score for switched models
pp_brier_excl2=scoring(pred=pp_preds_exclosure$preds_switch, response=pp_preds_exclosure$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine switched model predictions and brier score dataframes
pp_brier_excl_switch=cbind(pp_preds_exclosure, pp_brier_excl2)%>%rename(quadratic_switch=quadratic)

#calculate brier score differences
pp_brier_exclosure=left_join(pp_brier_excl_same, pp_brier_excl_switch)%>%mutate(treatment="removal",brier_diff=quadratic_same-quadratic_switch)

#select relevant rows with horizons 1,6 12
pp_brier_exclosure1=pp_brier_exclosure%>%group_by(m)%>%slice(which(row_number()%in%c(1,6,12)))

#subdivide into each horizon for easier plotting
ppbx1=pp_brier_exclosure1%>%group_by(m)%>%slice_head(n=1)%>%mutate(horizon="1")
ppbx6=pp_brier_exclosure1%>%group_by(m)%>%slice(which(row_number()==2))%>%mutate(horizon="6")
ppbx12=pp_brier_exclosure1%>%group_by(m)%>%slice(which(row_number()==3))%>%mutate(horizon="12")

#combine all 3 horizons
ppbx=rbind(ppbx1,ppbx6, ppbx12)

#combine control and exclosure and filter out per horizon for plotting
pp_briers1=rbind(ppb,ppbx)%>%filter(horizon==1)
pp_briers6=rbind(ppb,ppbx)%>%filter(horizon==6)
pp_briers12=rbind(ppb,ppbx)%>%filter(horizon==12)
