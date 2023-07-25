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
require(ggplot2)
require(ggpubr)
require(tidyverse)

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

#read covariates data####

covars=weather(level="newmoon", fill=TRUE, horizon=365, path=get_default_data_path())%>%
  select(newmoonnumber,meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)%>%
  mutate(meantemp_lag1=lag(meantemp,order_by=newmoonnumber))

ppcont_covs=right_join(covars,ppcont_dat)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PP)%>%rename("abundance"="PP")

ppexcl_covs=right_join(covars,ppexcl_dat)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PP)%>%rename("abundance"="PP")

#select data from 2010-2019
ppcont_covs_ro=ppcont_covs%>% filter(!newmoonnumber<403, !newmoonnumber>526)%>%
  mutate(part = ifelse(newmoonnumber<=476,"Train","Test")) #476 should not be hardcoded

ppexcl_covs_ro=ppexcl_covs%>% filter(!newmoonnumber<403, !newmoonnumber>526)%>%
  mutate(part = ifelse(newmoonnumber<=476,"Train","Test"))

ppcont_covs_ro$abundance=round_na.interp(ppcont_covs_ro$abundance)
ppexcl_covs_ro$abundance=round_na.interp(ppexcl_covs_ro$abundance)

#look at time-series (initial)
pp1=ggplot(data=ppcont_covs, aes(newmoonnumber, abundance)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("PP control abundances")

pp2=ggplot(data=ppexcl_covs, aes(newmoonnumber, abundance, color = part)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("PP exclosure abundances")

ggarrange(pp1, pp2, common.legend = T)

ggplot(data=ppcont_covs, aes(newmoonnumber, abundance)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("PP abundances")+annotate("rect", alpha=.2, fill='red',xmin=403, xmax=468, ymin=0, ymax=90)+
  annotate("rect", alpha=.2, fill='blue',xmin=469, xmax=526, ymin=0, ymax=90)

#create rolling origin object for analysis####
n_moons_yr=13
n_yrs=5
n_moons_train=n_moons_yr*n_yrs
n_moons_test=n_moons_yr*1

PPcontrol_dat <- 
  rolling_origin(
    data       = ppcont_covs_ro, #all PP control data (2010-2019)
    initial    = n_moons_train, #samples used for modelling (training)
    assess     = n_moons_test, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed; 
  ) 

PPexclosure_dat <- 
  rolling_origin(
    data       = ppexcl_covs_ro, #all PP exclosure data (2010-2019)
    initial    = n_moons_train, #samples used for modelling (training)
    assess     = n_moons_test, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

#add column for model(same data and model)####
PPcontrol_dat$model=map(PPcontrol_dat$splits, rolling_mod)
PPexclosure_dat$model=map(PPexclosure_dat$splits, rolling_mod)

#add column for model predictions (same data and model)####
PPcontrol_dat$preds_same=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$model, PPcontrol_dat$model), get_preds)
PPexclosure_dat$preds_same=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$model, PPexclosure_dat$model), get_preds)

#add column for model predictions from switched model####
PPcontrol_dat$preds_switch=pmap(list(PPcontrol_dat$splits, PPexclosure_dat$model, PPcontrol_dat$model), get_preds)
PPexclosure_dat$preds_switch=pmap(list(PPexclosure_dat$splits, PPcontrol_dat$model, PPexclosure_dat$model), get_preds)

#add column for model evals from same model####
PPcontrol_dat$evals_same=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$preds_same), mod_evals_same)
PPexclosure_dat$evals_same=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$preds_same), mod_evals_same)

#add column for model evals from switched model####
PPcontrol_dat$evals_switch=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$preds_switch),mod_evals_switch)
PPexclosure_dat$evals_switch=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$preds_switch),mod_evals_switch)

#interpolated data length####
mdate=read.csv("https://raw.githubusercontent.com/weecology/PortalData/main/Rodents/moon_dates.csv")

ppm=mdate%>%filter(!newmoonnumber<403, !newmoonnumber>526)
ppm_na=ppm%>%filter(is.na(censusdate))
length(ppm_na$newmoonnumber)/length(ppm$newmoonnumber)