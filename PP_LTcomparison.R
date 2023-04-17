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
require(ggplot2)
require(tidyverse)

#portalcasting setup####
main <- "~/portalcast_directory"

#read data####

pp_cont_interp=read_data(main = main, data_name = "rodents_table", dataset = "controls_interp")%>%
  select(newmoonnumber, PP)
pp_excl_interp=read_data(main = main, data_name = "rodents_table", dataset = "exclosures_interp")%>%
  select(newmoonnumber, PP)

#read covariates data####

covars=weather(level="newmoon", fill=TRUE, horizon=365, path=get_default_data_path())%>%
  select(newmoonnumber,meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)%>%
  mutate(meantemp_lag1=lag(meantemp,order_by=newmoonnumber))

ppcontrols_covs=right_join(covars,pp_cont_interp)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PP)%>%rename("abundance"="PP")

ppexcl_covs=right_join(covars,pp_excl_interp)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PP)%>%rename("abundance"="PP")

ppcont_dat=ppcontrols_covs%>% filter(!newmoonnumber<403, !newmoonnumber>526)%>%
  mutate(part = ifelse(newmoonnumber<=476,"Train","Test")) #476 should not be hardcoded

ppexcl_dat=ppexcl_covs%>% filter(!newmoonnumber<403, !newmoonnumber>526)%>%
  mutate(part = ifelse(newmoonnumber<=476,"Train","Test"))

#look at time-series (initial)
pp1=ggplot(data=ppcont_dat, aes(newmoonnumber, abundance, color = part)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("PP control abundances")

pp2=ggplot(data=ppexcl_dat, aes(newmoonnumber, abundance, color = part)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("PP exclosure abundances")

ggarrange(pp1, pp2, common.legend = T)

#create seasonality manually (CAN SKIP)
moons <- read_moons(main     = main,
                    settings = directory_settings())

moon_foys        <- foy(dates = as.Date(moons$newmoondate))
sin2pifoy        <- sin(2 * pi * moon_foys)
cos2pifoy        <- cos(2 * pi * moon_foys)
fouriers         <- data.frame(sin2pifoy, cos2pifoy)

match_time=which(moons$newmoonnumber %in% ppexcl_dat$newmoonnumber & moons$newmoonnumber)
fors=fouriers[match_time,]      

#create full data frame
pp_datc=cbind(fors, ppcont_dat)
pp_date=cbind(fors, ppexcl_dat)

#something for seasonal GARCH models?
past <- list(past_obs = c(1,13), external=TRUE) #autoregressive terms (1,13), external effect=T


#create rolling origin object for analysis####
n_moons_yr=13
n_yrs=5
n_moons_train=n_moons_yr*n_yrs
n_moons_test=n_moons_yr*1

PPcontrol_dat <- 
  rolling_origin(
    data       = pp_datc, #all PP control data (2010-2019)
    initial    = n_moons_train, #samples used for modelling (training)
    assess     = n_moons_test, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed; 
  ) 

PPexclosure_dat <- 
  rolling_origin(
    data       = pp_date, #all PP exclosure data (2010-2019)
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
