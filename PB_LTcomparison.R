#load packages####
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

#read data####
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

#read covariates data####

covars=weather(level="newmoon", fill=TRUE, horizon=365, path=get_default_data_path())%>%
  select(newmoonnumber,meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)%>%
  mutate(meantemp_lag1=lag(meantemp,order_by=newmoonnumber))

pbcont_covs=right_join(covars,pbcont_dat)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PB)%>%rename("abundance"="PB")

pbexcl_covs=right_join(covars,pbexcl_dat)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PB)%>%rename("abundance"="PB")

#select data from 2000-2010
pbcont_covs_ro=pbcont_covs%>%
  filter(!newmoonnumber<279, !newmoonnumber>414)%>%
  mutate(part = ifelse(newmoonnumber<=353,"Train","Test"))

pbexcl_covs_ro=pbexcl_covs%>%
  filter(!newmoonnumber<279, !newmoonnumber>414)%>%
  mutate(part = ifelse(newmoonnumber<=353,"Train","Test"))

pbcont_covs_ro$abundance=round_na.interp(pbcont_covs_ro$abundance)
pbexcl_covs_ro$abundance=round_na.interp(pbexcl_covs_ro$abundance)

#look at time-series
pb1=ggplot(data=pbcont_covs_ro, aes(newmoonnumber, abundance, color = part)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("PB control abundances")

pb2=ggplot(data=pbexcl_covs_ro, aes(newmoonnumber, abundance, color = part)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("PB exclosure abundances")

ggarrange(pb1, pb2, common.legend = T)

ggplot(data=pbcont_covs, aes(newmoonnumber, abundance)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("PB abundances")+annotate("rect", alpha=.2, fill='red',xmin=279, xmax=344, ymin=0, ymax=60)+
  annotate("rect", alpha=.2, fill='blue',xmin=345, xmax=414, ymin=0, ymax=60)

#something for seasonal GARCH models?
past <- list(past_obs = c(1,13), external=TRUE) #autoregressive terms (1,13), external effect=T

#create rolling origin object for analysis####
n_moons_yr=13
n_yrs=5
n_moons_train=n_moons_yr*n_yrs
n_moons_test=n_moons_yr*1

PBcontrol_dat <- 
  rolling_origin(
    data       = pbcont_covs_ro, #all PB control data (2000-2010)
    initial    = n_moons_train, #samples used for modelling (training)
    assess     = 13, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  ) 

PBexclosure_dat <- 
  rolling_origin(
    data       = pbexcl_covs_ro, #all PB exclosure data (2000-2010)
    initial    = n_moons_train, #samples used for modelling (training)
    assess     = 13, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

#add column for model(same data and model)####

PBcontrol_dat$model=map(PBcontrol_dat$splits, rolling_mod)
PBexclosure_dat$model=map(PBexclosure_dat$splits, rolling_mod)

#add column for model predictions (same data and model)####
PBcontrol_dat$preds_same=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$model, PBcontrol_dat$model), get_preds)
PBexclosure_dat$preds_same=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$model, PBexclosure_dat$model), get_preds)

#add column for model predictions from switched model####
PBcontrol_dat$preds_switch=pmap(list(PBcontrol_dat$splits, PBexclosure_dat$model, PBcontrol_dat$model), get_preds)
PBexclosure_dat$preds_switch=pmap(list(PBexclosure_dat$splits, PBcontrol_dat$model, PBexclosure_dat$model), get_preds)

#add column for model evals from same model####
PBcontrol_dat$evals_same=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$preds_same), mod_evals_same)
PBexclosure_dat$evals_same=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$preds_same), mod_evals_same)

#add column for model evals from switched model (h=1)####
PBcontrol_dat$evals_switch=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$preds_switch),mod_evals_switch)
PBexclosure_dat$evals_switch=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$preds_switch),mod_evals_switch)
