#script for conducting model transfer analyses###

#function to run nbevGARCH models
rolling_mod=function(split) {
  
  analysis_set= analysis(split) #get dataframe
  
  fit_model= tsglm(analysis_set[,"abundance"], model = list(past_obs=c(1,12), external=TRUE), 
                   distr = "nbinom", 
                   xreg  = analysis_set[,3:5], 
                   link  = "log")
}

#function to generate predictions for each model split
get_preds=function(split, model,y_source) {
  
  analysis_set= analysis(split) #get training data
  assessment_set=assessment(split) #get testing data
  model$ts = y_source$ts #select the initial condition from matching model
  model$response = y_source$response
  preds=predict(model, n.ahead=12, newxreg= assessment_set[,3:5])$pred
}

#function for evaluating model
mod_evals_same=function(split, preds_same) {
  
  assessment_set= assessment(split) #get validation data
  
  holdout= assessment_set[,"abundance"] #select relevant column in dataframe
  
  preds_same= preds_same
  
  rmse1_same=pmap(list(holdout[1], preds_same[1]), Metrics::rmse)
  rmse6_same=pmap(list(holdout[6], preds_same[6]), Metrics::rmse)
  rmse12_same=pmap(list(holdout[12], preds_same[12]), Metrics::rmse)
  
  rmse_same=cbind(rmse1_same, rmse6_same, rmse12_same)
}

mod_evals_switch=function(split, preds_switch) {
  
  assessment_set= assessment(split) #get validation data
  
  holdout= assessment_set[,"abundance"] #select relevant column in dataframe
  
  preds_switch= preds_switch
  
  rmse1_switch=pmap(list(holdout[1], preds_switch[1]), Metrics::rmse)
  rmse6_switch=pmap(list(holdout[6], preds_switch[6]), Metrics::rmse)
  rmse12_switch=pmap(list(holdout[12], preds_switch[12]), Metrics::rmse)
  
  rmse_switch=cbind(rmse1_switch, rmse6_switch, rmse12_switch)
}


#function for plotting forecasts
get_dat_same=function(split, preds_same) {
  
  assessment_set= assessment(split) #get validation data
  
  holdout= assessment_set[,"abundance"]
  moon= assessment_set[,"newmoonnumber"]
  preds_same= preds_same
  
  hp=cbind(moon, holdout, preds_same)
}

get_dat_switch=function(split, preds_switch) {
  
  assessment_set= assessment(split) #get validation data
  
  holdout= assessment_set[,"abundance"]
  moon= assessment_set[,"newmoonnumber"]
  preds_switch= preds_switch
  
  hp=cbind(moon, holdout, preds_switch)
}

#function for plotting forecast evals
get_evals1_diff=function(split, evals_same, evals_switch, id) {
  
  assessment_set= assessment(split) #get validation data
  
  newmoon= assessment_set[,"newmoonnumber"][1]
  
  score_same=as.vector(as.numeric(evals_same[1]))
  score_switch=as.vector(as.numeric(evals_switch[1]))
  score_diff=score_same-score_switch
  
  id=id
  h=1
  hp=cbind(newmoon,id, h, score_same, score_switch, score_diff)
}

get_evals6_diff=function(split, evals_same, evals_switch, id) {
  
  assessment_set= assessment(split) #get validation data
  
  newmoon= assessment_set[,"newmoonnumber"][6]
  
  score_same=as.vector(as.numeric(evals_same[2]))
  score_switch=as.vector(as.numeric(evals_switch[2]))
  score_diff=score_same-score_switch
  
  id=id
  h=6
  hp=cbind(newmoon,id, h, score_same, score_switch, score_diff)
}

get_evals12_diff=function(split, evals_same, evals_switch, id) {
  
  assessment_set= assessment(split) #get validation data
  
  newmoon= assessment_set[,"newmoonnumber"][12]
  
  score_same=as.vector(as.numeric(evals_same[3]))
  score_switch=as.vector(as.numeric(evals_switch[3]))
  score_diff=score_same-score_switch
  
  id=id
  h=12
  hp=cbind(newmoon,id, h, score_same, score_switch, score_diff)
}
