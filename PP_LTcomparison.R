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

#load Portal data####
portal_dat1=summarise_rodent_data(
  path = get_default_data_path(),
  clean = TRUE,
  level = "Treatment",
  type = "Rodents",
  unknowns = FALSE,
  shape = "flat",
  time = "all",
  output="abundance"
)%>%mutate(year=lubridate::year(censusdate),
           month=lubridate::month(censusdate),
           day=lubridate::day(censusdate))%>%
  filter(!(treatment%in%c("removal","spectabs")))%>%
  relocate(c("year", "month", "day"), .before="treatment")

#get covariates
temp=weather(level="daily", fill=TRUE, horizon=90)

#temp$date=as.Date(paste(temp$year, temp$month, 01), "%Y %m %d")

#manipulate data####
#create dataframe with counts and covariates
count_covs=right_join(temp,portal_dat1, by=c("year", "month", "day"))%>%
  select(date,censusdate, year, month, day, treatment, mintemp, maxtemp, meantemp, precipitation,
         warm_precip, cool_precip, period, species, abundance)

#create distinction between training and test data
count_dat=count_covs%>% filter(!year<1995, species=="PP")%>%
  mutate(part = ifelse(year <= 2015,"Train","Test"))

ggplot(data=count_dat, aes(period, abundance, color = part)) +
  geom_point(alpha = 0.2) +
  geom_line()

#create seasonality
moon_foys        <- foy(dates = as.Date(count_dat$censusdate))
sin2pifoy        <- sin(2 * pi * moon_foys)
cos2pifoy        <- cos(2 * pi * moon_foys)
fouriers         <- data.frame(sin2pifoy, cos2pifoy, count_dat$meantemp, count_dat$precipitation)

#select portion of dataframe for control and exclosure and for predictions

port_moonc=which(count_dat$treatment=="control" & count_dat$part=="Train")
port_moone=which(count_dat$treatment=="exclosure"& count_dat$part=="Train")
test_dat1=which(count_dat$part=="Test" & count_dat$treatment=="control")
test_dat2=which(count_dat$part=="Test" & count_dat$treatment=="exclosure")

predsev1=fouriers[port_moonc,]
predsev2=fouriers[port_moone,]
predsev3=fouriers[test_dat1,]
predsev3=fouriers[test_dat2,]

#something for seasonal GARCH models?
past <- list(past_obs = 1, past_mean = 12) #p = 1 (first-order AR),q = 12 (approximately yearly moving average)

#get all pp control data
pp_datc=count_dat%>%filter(treatment=="control")

#get all pp exclosure data (1995-2020)
pp_date=count_dat%>%filter(treatment=="exclosure")

#get control training data (1995-2015)
ppm_datc=count_dat[port_moonc,]

#get exclosure training ata
ppm_date=count_dat[port_moone,]

#get control test data (2016-2019)
ppm_datc2=count_dat[test_dat1,]

#get exclosure training ata
ppm_date2=count_dat[test_dat2,]

#create dataframe for analysis####
PPcontrol_dat <- 
  rolling_origin(
    data       = pp_datc, #all PP control data (1995-2020)
    initial    = 120, #samples used for modelling (training)
    assess     = 12, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  ) 

PPexclosure_dat <- 
  rolling_origin(
    data       = pp_date, #all PP control data (1995-2020)
    initial    = 120, #samples used for modelling (training)
    assess     = 12, # number of samples used for each assessment resample (horizon)
    cumulative = FALSE #length of analysis set is fixed
  )

#check what the dataframe looks like
#t1=PPcontrol_dat$splits[[2]] 
#t2=PPexclosure_dat$splits[[2]
#analysis(t1) %>% 
#  tail()

#assessment(t1) %>% 
#  tail()

source("D:/WeecologyProjects/portalcasting/Ltcompare_functions.R")

#add column for model(same data and model)####
PPcontrol_dat$model=map(PPcontrol_dat$splits, rolling_mod)
PPexclosure_dat$model=map(PPexclosure_dat$splits, rolling_mod)

#add column for model coefs, CAN SKIP THIS####
#PPcontrol_dat$modcoef_control=map(PPcontrol_dat$splits, rolling_mod_coef)
#PPexclosure_dat$modcoef_exclosure=map(PPexclosure_dat$splits, rolling_mod_coef)

#add column for model predictions (same data and model)####
PPcontrol_dat$preds_same=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$model), get_preds)
PPexclosure_dat$preds_same=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$model), get_preds)

#add column for model predictions from switched model####
PPcontrol_dat$preds_switch=pmap(list(PPcontrol_dat$splits, PPexclosure_dat$model), get_preds)
PPexclosure_dat$preds_switch=pmap(list(PPexclosure_dat$splits, PPcontrol_dat$model), get_preds)

#add column for model evals from same model(h=1)####
PPcontrol_dat$evals_same=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$model), mod_evals1)
PPexclosure_dat$evals_same=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$model), mod_evals1)

#add column for model evals from switched model (h=1)####
PPcontrol_dat$evals_switch=pmap(list(PPcontrol_dat$splits,PPexclosure_dat$model),mod_evals1)
PPexclosure_dat$evals_switch=pmap(list(PPexclosure_dat$splits,PPcontrol_dat$model),mod_evals1)

#add column for model evals from same model(h=6)####
PPcontrol_dat$evals_same6=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$model), mod_evals6)
PPexclosure_dat$evals_same6=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$model), mod_evals6)

#add column for model evals from switched model (h=6)####
PPcontrol_dat$evals_switch6=pmap(list(PPcontrol_dat$splits,PPexclosure_dat$model),mod_evals6)
PPexclosure_dat$evals_switch6=pmap(list(PPexclosure_dat$splits,PPcontrol_dat$model),mod_evals6)

#add column for model evals from same model(h=12)####
PPcontrol_dat$evals_same12=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$model), mod_evals12)
PPexclosure_dat$evals_same12=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$model), mod_evals12)

#add column for model evals from switched model (h=12)####
PPcontrol_dat$evals_switch12=pmap(list(PPcontrol_dat$splits,PPexclosure_dat$model),mod_evals12)
PPexclosure_dat$evals_switch12=pmap(list(PPexclosure_dat$splits,PPcontrol_dat$model),mod_evals12)

