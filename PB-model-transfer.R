#Script to run model transfer analyses: PB case study###

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

pbcont_dat=rodent_data%>%
  select(newmoonnumber,treatment, PB)%>%
  replace_na(list(treatment='control'))%>%
  filter(treatment=="control")

pbexcl_dat=rodent_data%>%
  select(newmoonnumber,treatment, PB)%>%
  replace_na(list(treatment='exclosure'))%>%
  filter(treatment=="exclosure")

#covariates data

covars=weather(level="newmoon", fill=TRUE, horizon=365, path=get_default_data_path())%>%
  select(newmoonnumber,meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)%>%
  mutate(meantemp_lag1=lag(meantemp,order_by=newmoonnumber))

pbcont_covs=right_join(covars,pbcont_dat)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PB)%>%rename("abundance"="PB")

pbexcl_covs=right_join(covars,pbexcl_dat)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PB)%>%rename("abundance"="PB")

#select data from Dec 1999- June 2009
pbcont_covs_dat=pbcont_covs%>%filter(!newmoonnumber<278, !newmoonnumber>396)
pbexcl_covs_dat=pbexcl_covs%>%filter(!newmoonnumber<278, !newmoonnumber>396)

#interpolate missing data
pbcont_covs_dat$abundance=round_na.interp(pbcont_covs_dat$abundance)
pbexcl_covs_dat$abundance=round_na.interp(pbexcl_covs_dat$abundance)

#create rolling origin object for analysis
n_moons_yr=12
n_yrs=5
n_moons_train=n_moons_yr*n_yrs
n_moons_test=n_moons_yr*1

PBcontrol_dat <- 
  rolling_origin(
    data       = pbcont_covs_dat, #all PB control data (1999-2009)
    initial    = n_moons_train, #samples used for modelling (training)
    assess     = n_moons_yr, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  ) 

PBexclosure_dat <- 
  rolling_origin(
    data       = pbexcl_covs_dat, #all PB exclosure data (1999-2009)
    initial    = n_moons_train, #samples used for modelling (training)
    assess     = n_moons_yr, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

#DATA ANALYSES####

#Model fitting####
#fitting matched and mismatched models

#add column for model(same data and model)
PBcontrol_dat$model=map(PBcontrol_dat$splits, rolling_mod)
PBexclosure_dat$model=map(PBexclosure_dat$splits, rolling_mod)

#matched and mismatched model predictions

#add column for model predictions (same data and model)
PBcontrol_dat$preds_same=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$model, PBcontrol_dat$model), get_preds)
PBexclosure_dat$preds_same=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$model, PBexclosure_dat$model), get_preds)

#add column for model predictions from switched model
PBcontrol_dat$preds_switch=pmap(list(PBcontrol_dat$splits, PBexclosure_dat$model, PBcontrol_dat$model), get_preds)
PBexclosure_dat$preds_switch=pmap(list(PBexclosure_dat$splits, PBcontrol_dat$model, PBexclosure_dat$model), get_preds)

#matched and mismatched model evaluations

#add column for model evals from same model
PBcontrol_dat$evals_same=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$preds_same), mod_evals_same)
PBexclosure_dat$evals_same=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$preds_same), mod_evals_same)

#add column for model evals from switched model
PBcontrol_dat$evals_switch=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$preds_switch),mod_evals_switch)
PBexclosure_dat$evals_switch=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$preds_switch),mod_evals_switch)

#Parameter comparison####

pbcontrol_int_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(1)
pbexclosure_int_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(1)

pbcontrol_b1_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(2)
pbexclosure_b1_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(2)

pbcontrol_b12_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(3)
pbexclosure_b12_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(3)

pbcontrol_temp1_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(4)
pbexclosure_temp1_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(4)

pbcontrol_warmprec_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(5)
pbexclosure_warmprec_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(5)

pbcontrol_coolprec_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(6)
pbexclosure_coolprec_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(6)

pbwarmprec=cbind(pbcontrol_warmprec_coefs, pbexclosure_warmprec_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_warmprec_coefs", "removal"="pbexclosure_warmprec_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "warm_precip")

pbcoolprec=cbind(pbcontrol_coolprec_coefs, pbexclosure_coolprec_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_coolprec_coefs", "removal"="pbexclosure_coolprec_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "cool_precip")

pbtemps=cbind(pbcontrol_temp1_coefs, pbexclosure_temp1_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_temp1_coefs", "removal"="pbexclosure_temp1_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "temp")

pbints=cbind(pbcontrol_int_coefs, pbexclosure_int_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_int_coefs", "removal"="pbexclosure_int_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "intercept")

pbb1=cbind(pbcontrol_b1_coefs, pbexclosure_b1_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_b1_coefs", "removal"="pbexclosure_b1_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "beta1")

pbb12=cbind(pbcontrol_b12_coefs, pbexclosure_b12_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_b12_coefs", "removal"="pbexclosure_b12_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "beta12")

coef_df_PB=as.data.frame(list(pbints, pbb1, pbb12, pbtemps, pbwarmprec, pbcoolprec))%>%
  select(treatment, intercept, beta1, beta12, temp,cool_precip, warm_precip)

int_cont=coef_df_PB%>%filter(treatment=="control")%>%select(intercept)%>%rename("control"="intercept")
int_excl=coef_df_PB%>%filter(treatment=="removal")%>%select(intercept)%>%rename("removal"="intercept")

ar1_cont=coef_df_PB%>%filter(treatment=="control")%>%select(beta1)%>%rename("control"="beta1")
ar1_excl=coef_df_PB%>%filter(treatment=="removal")%>%select(beta1)%>%rename("removal"="beta1")

ar12_cont=coef_df_PB%>%filter(treatment=="control")%>%select(beta12)%>%rename("control"="beta12")
ar12_excl=coef_df_PB%>%filter(treatment=="removal")%>%select(beta12)%>%rename("removal"="beta12")

temp_cont=coef_df_PB%>%filter(treatment=="control")%>%select(temp)%>%rename("control"="temp")
temp_excl=coef_df_PB%>%filter(treatment=="removal")%>%select(temp)%>%rename("removal"="temp")

cprec_cont=coef_df_PB%>%filter(treatment=="control")%>%select(cool_precip)%>%rename("control"="cool_precip")
cprec_excl=coef_df_PB%>%filter(treatment=="removal")%>%select(cool_precip)%>%rename("removal"="cool_precip")

wprec_cont=coef_df_PB%>%filter(treatment=="control")%>%select(warm_precip)%>%rename("control"="warm_precip")
wprec_excl=coef_df_PB%>%filter(treatment=="removal")%>%select(warm_precip)%>%rename("removal"="warm_precip")

o1=cbind(int_cont, int_excl)
o2=cbind(ar1_cont, ar1_excl)
o3=cbind(ar12_cont, ar12_excl)
o4=cbind(temp_cont, temp_excl)
o5=cbind(cprec_cont, cprec_excl)
o6=cbind(wprec_cont, wprec_excl)

###directional shift####
length(which(o1>0))/96
length(which(o2>0))/96
length(which(o3>0))/96
length(which(o4>0))/96
length(which(o5>0))/96
length(which(o6>0))/96

##degree of overlap####
overlap(o1, plot=T)
overlap(o2, plot=T)
overlap(o3, plot=T)
overlap(o4, plot=T)
overlap(o5, plot=T)
overlap(o6, plot=T)

##multiple comparisons test####

#calculate pvalue following: https://www.bmj.com/content/343/bmj.d2304 

PBcontrol_dat$pvalue=map(PBcontrol_dat$model, get_pvalue)
PBexclosure_dat$pvalue=map(PBexclosure_dat$model, get_pvalue)

#select pvalues for each parameter

pbcont_int_pval=PBcontrol_dat$pvalue%>%map(unlist)%>%map_dbl(1)
pbcont_b1_pval=PBcontrol_dat$pvalue%>%map(unlist)%>%map_dbl(2)
pbcont_b12_pval=PBcontrol_dat$pvalue%>%map(unlist)%>%map_dbl(3)
pbcont_temp_pval=PBcontrol_dat$pvalue%>%map(unlist)%>%map_dbl(4)
pbcont_wprec_pval=PBcontrol_dat$pvalue%>%map(unlist)%>%map_dbl(5)
pbcont_cprec_pval=PBcontrol_dat$pvalue%>%map(unlist)%>%map_dbl(6)

pbexcl_int_pval=PBexclosure_dat$pvalue%>%map(unlist)%>%map_dbl(1)
pbexcl_b1_pval=PBexclosure_dat$pvalue%>%map(unlist)%>%map_dbl(2)
pbexcl_b12_pval=PBexclosure_dat$pvalue%>%map(unlist)%>%map_dbl(3)
pbexcl_temp_pval=PBexclosure_dat$pvalue%>%map(unlist)%>%map_dbl(4)
pbexcl_wprec_pval=PBexclosure_dat$pvalue%>%map(unlist)%>%map_dbl(5)
pbexcl_cprec_pval=PBexclosure_dat$pvalue%>%map(unlist)%>%map_dbl(6)

#combine dataframes

pbint_pval=cbind(pbcont_int_pval, pbexcl_int_pval)%>%as.data.frame%>%
  rename("control"="pbcont_int_pval", "removal"="pbexcl_int_pval")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "raw_pvalue")

pbb1_pval=cbind(pbcont_b1_pval, pbexcl_b1_pval)%>%as.data.frame%>%
  rename("control"="pbcont_b1_pval", "removal"="pbexcl_b1_pval")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "raw_pvalue")

pbb12_pval=cbind(pbcont_b12_pval, pbexcl_b12_pval)%>%as.data.frame%>%
  rename("control"="pbcont_b12_pval", "removal"="pbexcl_b12_pval")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "raw_pvalue")

pbtemp_pval=cbind(pbcont_temp_pval, pbexcl_temp_pval)%>%as.data.frame%>%
  rename("control"="pbcont_temp_pval", "removal"="pbexcl_temp_pval")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "raw_pvalue")

pbwarmprec_pval=cbind(pbcont_wprec_pval, pbexcl_wprec_pval)%>%as.data.frame%>%
  rename("control"="pbcont_wprec_pval", "removal"="pbexcl_wprec_pval")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "raw_pvalue")

pbcoolprec_pval=cbind(pbcont_cprec_pval, pbexcl_cprec_pval)%>%as.data.frame%>%
  rename("control"="pbcont_cprec_pval", "removal"="pbexcl_cprec_pval")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "raw_pvalue")

#pairwise t-test
pairwise.t.test(pbint_pval$raw_pvalue, pbint_pval$treatment, p.adjust.method="fdr", pool.sd=FALSE, paired=TRUE) 
pairwise.t.test(pbb1_pval$raw_pvalue, pbb1_pval$treatment, p.adjust.method="fdr", pool.sd=FALSE, paired=TRUE) 
pairwise.t.test(pbb12_pval$raw_pvalue, pbb12_pval$treatment, p.adjust.method="fdr", pool.sd=FALSE, paired=TRUE) 
pairwise.t.test(pbtemp_pval$raw_pvalue, pbtemp_pval$treatment, p.adjust.method="fdr", pool.sd=FALSE, paired=TRUE) 
pairwise.t.test(pbwarmprec_pval$raw_pvalue, pbwarmprec_pval$treatment, p.adjust.method="fdr",pool.sd=FALSE, paired=TRUE) 
pairwise.t.test(pbcoolprec_pval$raw_pvalue, pbcoolprec_pval$treatment, p.adjust.method="fdr", pool.sd=FALSE, paired=TRUE) 

#add adjusted p-values (FDR)
pbint_pval$FDR=p.adjust(pbint_pval$raw_pvalue, method="fdr")
pbb1_pval$FDR=p.adjust(pbb1_pval$raw_pvalue, method="fdr")
pbb12_pval$FDR=p.adjust(pbb12_pval$raw_pvalue, method="fdr")
pbtemp_pval$FDR=p.adjust(pbtemp_pval$raw_pvalue, method="fdr")
pbwarmprec_pval$FDR=p.adjust(pbwarmprec_pval$raw_pvalue, method="fdr")
pbcoolprec_pval$FDR=p.adjust(pbcoolprec_pval$raw_pvalue, method="fdr")

#Evaluate model transferability####

###generate forecasts####

#control-control
pbcont_preds_same=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$preds_same), get_dat_same)

#control dat-exclosure mod
pbcont_preds_switch=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$preds_switch), get_dat_switch)

#exclosure-exclosure
pbexcl_preds_same=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$preds_same), get_dat_same)

#exclosure dat-control mod
pbexcl_preds_switch=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$preds_switch), get_dat_switch)

m=rep(seq(1:48), each=12)
code1="same"
code2="switched"

PBpreds_cont_same=do.call(rbind.data.frame, pbcont_preds_same)
PBpreds_cont_same=cbind(PBpreds_cont_same, m, code1)
PBpreds_cont_switch=do.call(rbind.data.frame, pbcont_preds_switch)
PBpreds_cont_switch=cbind(PBpreds_cont_switch, m, code2)

PBpreds_excl_same=do.call(rbind.data.frame, pbexcl_preds_same)
PBpreds_excl_same=cbind(PBpreds_excl_same, m, code1)
PBpreds_excl_switch=do.call(rbind.data.frame, pbexcl_preds_switch)
PBpreds_excl_switch=cbind(PBpreds_excl_switch, m, code2)

###RMSE####

#h=1
pbcont_evals1_diff=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$evals_same, PBcontrol_dat$evals_switch ,PBcontrol_dat$id), get_evals1_diff)
pbevals1_cont_diff=do.call(rbind.data.frame, pbcont_evals1_diff)%>%mutate(plot="control")

pbexcl_evals1_diff=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$evals_same, PBexclosure_dat$evals_switch ,PBexclosure_dat$id), get_evals1_diff)
pbevals1_excl_diff=do.call(rbind.data.frame, pbexcl_evals1_diff)%>%mutate(plot="removal")

pbevals1=rbind(pbevals1_cont_diff, pbevals1_excl_diff)
pbevals1$newmoon=as.integer(pbevals1$newmoon)
pbevals1$h=as.integer(pbevals1$h)
pbevals1$score_same=as.numeric(pbevals1$score_same)
pbevals1$score_switch=as.numeric(pbevals1$score_switch)
pbevals1$score_diff=as.numeric(pbevals1$score_diff)

#h=6
pbcont_evals6_diff=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$evals_same, PBcontrol_dat$evals_switch ,PBcontrol_dat$id), get_evals6_diff)
pbevals6_cont_diff=do.call(rbind.data.frame, pbcont_evals6_diff)%>%mutate(plot="control")

pbexcl_evals6_diff=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$evals_same, PBexclosure_dat$evals_switch ,PBexclosure_dat$id), get_evals6_diff)
pbevals6_excl_diff=do.call(rbind.data.frame, pbexcl_evals6_diff)%>%mutate(plot="removal")

pbevals6=rbind(pbevals6_cont_diff, pbevals6_excl_diff)
pbevals6$newmoon=as.integer(pbevals6$newmoon)
pbevals6$h=as.integer(pbevals6$h)
pbevals6$score_same=as.numeric(pbevals6$score_same)
pbevals6$score_switch=as.numeric(pbevals6$score_switch)
pbevals6$score_diff=as.numeric(pbevals6$score_diff)

#h=12
pbcont_evals12_diff=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$evals_same, PBcontrol_dat$evals_switch ,PBcontrol_dat$id), get_evals12_diff)
pbevals12_cont_diff=do.call(rbind.data.frame, pbcont_evals12_diff)%>%mutate(plot="control")

pbexcl_evals12_diff=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$evals_same, PBexclosure_dat$evals_switch ,PBexclosure_dat$id), get_evals12_diff)
pbevals12_excl_diff=do.call(rbind.data.frame, pbexcl_evals12_diff)%>%mutate(plot="removal")

pbevals12=rbind(pbevals12_cont_diff, pbevals12_excl_diff)
pbevals12$newmoon=as.integer(pbevals12$newmoon)
pbevals12$h=as.integer(pbevals12$h)
pbevals12$score_same=as.numeric(pbevals12$score_same)
pbevals12$score_switch=as.numeric(pbevals12$score_switch)
pbevals12$score_diff=as.numeric(pbevals12$score_diff)

##Brier scores####

#control###

#combine predictions on same and switched models for control data
pb_preds_control=left_join(PBpreds_cont_same, PBpreds_cont_switch, by=c("moon", "holdout", "m"))

#calculate brier score for same models
pb_brier_cont1=scoring(pred=pb_preds_control$preds_same, response=pb_preds_control$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine same model predictions and brier score dataframes
pb_brier_cont_same=cbind(pb_preds_control, pb_brier_cont1)%>%rename(quadratic_same=quadratic)

#calculate brier score for switched models
pb_brier_cont2=scoring(pred=pb_preds_control$preds_switch, response=pb_preds_control$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine switched model predictions and brier score dataframes
pb_brier_cont_switch=cbind(pb_preds_control, pb_brier_cont2)%>%rename(quadratic_switch=quadratic)

#calculate brier score differences
pb_brier_control=left_join(pb_brier_cont_same, pb_brier_cont_switch)%>%mutate(treatment="control",brier_diff=quadratic_same-quadratic_switch)

#select relevant rows with horizons 1,6 12
pb_brier_control1=pb_brier_control%>%group_by(m)%>%slice(which(row_number()%in%c(1,6,12)))

#subdivide into each horizon for easier plotting
pbh1=pb_brier_control1%>%group_by(m)%>%slice_head(n=1)
pbh6=pb_brier_control1%>%group_by(m)%>%slice(which(row_number()==2))%>%mutate(horizon="6")
pbh12=pb_brier_control1%>%group_by(m)%>%slice(which(row_number()==3))%>%mutate(horizon="12")

#combine all 3 horizons
pbh=rbind(pbh1,pbh6,pbh12)

#exclosure###

#combine predictions on same and switched models for exclosure data
pb_preds_exclosure=left_join(PBpreds_excl_same, PBpreds_excl_switch, by=c("moon", "holdout", "m"))

#calculate brier score for same models
pb_brier_excl1=scoring(pred=pb_preds_exclosure$preds_same, response=pb_preds_exclosure$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine same model predictions and brier score dataframes
pb_brier_excl_same=cbind(pb_preds_exclosure, pb_brier_excl1)%>%rename(quadratic_same=quadratic)

#calculate brier score for switched models
pb_brier_excl2=scoring(pred=pb_preds_exclosure$preds_switch, response=pb_preds_exclosure$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine switched model predictions and brier score dataframes
pb_brier_excl_switch=cbind(pb_preds_exclosure, pb_brier_excl2)%>%rename(quadratic_switch=quadratic)

#calculate brier score differences
pb_brier_exclosure=left_join(pb_brier_excl_same, pb_brier_excl_switch)%>%mutate(treatment="removal",brier_diff=quadratic_same-quadratic_switch)

#select relevant rows with horizons 1,6 12
pb_brier_exclosure1=pb_brier_exclosure%>%group_by(m)%>%slice(which(row_number()%in%c(1,6,12)))

#subdivide into each horizon for easier plotting
pbbx1=pb_brier_exclosure1%>%group_by(m)%>%slice_head(n=1)%>%mutate(horizon="1")
pbbx6=pb_brier_exclosure1%>%group_by(m)%>%slice(which(row_number()==2))%>%mutate(horizon="6")
pbbx12=pb_brier_exclosure1%>%group_by(m)%>%slice(which(row_number()==3))%>%mutate(horizon="12")

#combine all 3 horizons
pbbx=rbind(pbbx1,pbbx6, pbbx12)

#combine control and exclosure and filter out per horizon for plotting
pb_briers1=rbind(pbb,pbbx)%>%filter(horizon==1)
pb_briers6=rbind(pbb,pbbx)%>%filter(horizon==6)
pb_briers12=rbind(pbb,pbbx)%>%filter(horizon==12)
