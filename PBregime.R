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

#select data from Dec 1999- Dec 2010
pbcont_covs_dat=pbcont_covs%>%filter(!newmoonnumber<278, !newmoonnumber>414)

#interpolate missing data
pbcont_covs_dat$abundance=round_na.interp(pbcont_covs_dat$abundance)

#rolling origin object for analysis####
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

#DATA ANALYSES####

#fitting matched and mismatched models

#add column for model(same data and model)
PBcontrol_dat$model=map(PBcontrol_dat$splits, rolling_mod)

#add column for model predictions (same data and model)
PBcontrol_dat$preds_same=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$model, PBcontrol_dat$model), get_preds)

#generate forecasts
#control-control
pbcont_preds_same=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$preds_same), get_dat_same)


#DATA VIZ####
pbpreds_plot1=ggplot()+
  geom_line(data=pbcont_covs_dat, aes(y=abundance, x=newmoonnumber), size=0.75)+
  geom_line(data=PBpreds_cont_same, aes(y=preds_same, x=moon, group=m, col="red"), alpha=0.4)+
  xlim(330,415)+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_rect(fill="grey",alpha=0.5,
            aes(ymin=0, ymax=700, xmin=409, xmax=414))+
  ggtitle("control")+xlab("sampling number (new moon number)")
 
pbpreds_plot1
annotate_figure(pbpreds_plot1,left = text_grob("C. baileyi", face = "italic", size = 16, rot=90))
