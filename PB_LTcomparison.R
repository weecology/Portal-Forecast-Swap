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

pb_cont_interp=read_data(main = main, data_name = "rodents_table", dataset = "controls_interp")%>%
  select(newmoonnumber, PB)
pb_excl_interp=read_data(main = main, data_name = "rodents_table", dataset = "exclosures_interp")%>%
  select(newmoonnumber, PB)

#read covariates data####

covars=weather(level="newmoon", fill=TRUE, horizon=365, path=get_default_data_path())%>%
  select(newmoonnumber,meantemp, mintemp, maxtemp, precipitation, warm_precip, cool_precip)%>%
  mutate(meantemp_lag1=lag(meantemp,order_by=newmoonnumber))

pbcontrols_covs=right_join(covars,pb_cont_interp)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PB)%>%rename("abundance"="PB")

pbexcl_covs=right_join(covars,pb_excl_interp)%>%
  select(newmoonnumber, meantemp, meantemp_lag1,
         warm_precip, cool_precip, PB)%>%rename("abundance"="PB")

#select data from 2000-2010
pbcont_dat=pbcontrols_covs%>% filter(!newmoonnumber<279, !newmoonnumber>414)%>%
  mutate(part = ifelse(newmoonnumber<=353,"Train","Test"))

pbexcl_dat=pbexcl_covs%>% filter(!newmoonnumber<279, !newmoonnumber>414)%>%
  mutate(part = ifelse(newmoonnumber<=353,"Train","Test"))

#look at time-series
pb1=ggplot(data=pbcont_dat, aes(newmoonnumber, abundance, color = part)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("PB control abundances")

pb2=ggplot(data=pbexcl_dat, aes(newmoonnumber, abundance, color = part)) +
  geom_point(alpha = 0.5, pch=19) +theme_classic()+geom_line()+
  ggtitle("PB exclosure abundances")

ggarrange(pb1, pb2, common.legend = T)

#create seasonality
moons <- read_moons(main     = main,
                    settings = directory_settings())

moon_foys        <- foy(dates = as.Date(moons$newmoondate))
sin2pifoy        <- sin(2 * pi * moon_foys)
cos2pifoy        <- cos(2 * pi * moon_foys)
fouriers         <- data.frame(sin2pifoy, cos2pifoy)

match_time=which(moons$newmoonnumber %in% pbexcl_dat$newmoonnumber & moons$newmoonnumber)
fors=fouriers[match_time,]      

#create full data frame
pb_datc=cbind(fors, pbcont_dat)
pb_date=cbind(fors, pbexcl_dat)

#something for seasonal GARCH models?
past <- list(past_obs = c(1,13), external=TRUE) #autoregressive terms (1,13), external effect=T

#create rolling origin object for analysis####
n_moons_yr=13
n_yrs=5
n_moons_train=n_moons_yr*n_yrs
n_moons_test=n_moons_yr*1

PBcontrol_dat <- 
  rolling_origin(
    data       = pb_datc, #all PB control data (2000-2010)
    initial    = n_moons_train, #samples used for modelling (training)
    assess     = 13, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  ) 

PBexclosure_dat <- 
  rolling_origin(
    data       = pb_date, #all PB exclosure data (2000-2010)
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
