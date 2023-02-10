#LTcompare_functions

#create function to run nbevGARCH models
rolling_mod=function(split) {
  
  analysis_set= analysis(split) #get dataframe
  
  fit_model= tsglm(analysis_set[,"abundance"], model = past, distr = "nbinom", 
                   xreg  = analysis_set[,5:7], 
                   link  = "log")
}


#create function to get model coefficients

rolling_mod_coef= function (split) {
  
  analysis_set= analysis(split) #get dataframe
  
  fit_model= tsglm(analysis_set[,"abundance"], model = past, distr = "nbinom", 
                   xreg  = analysis_set[,5:7], 
                   link  = "log")
  
  #get model coefficients    
  mod_coef=fit_model%>%
    coef()%>%
    tidy()
  
}

#create function to generate predictions for each model split
get_preds=function(split, model) {
  
  analysis_set= analysis(split) #get dataframe
  assessment_set=assessment(split)
  
  preds=predict(model, n.ahead=12, newxreg= assessment_set[,5:7])$pred
}

#create function for evaluating model
mod_evals1=function(split, model) {
  
  assessment_set= assessment(split) #get validation data
  
  holdout= assessment_set[,"abundance"] #select relevant column in dataframe
  
  mod_preds= predict(model, n.ahead=1)$pred
  
  rmse=pmap(list(holdout[1], mod_preds[1]), Metrics::rmse)
  
}

mod_evals6=function(split, model) {
  
  assessment_set= assessment(split) #get validation data
  
  holdout= assessment_set[,"abundance"] #select relevant column in dataframe
  
  mod_preds= predict(model, n.ahead=6)$pred
  
  rmse=pmap(list(holdout[6], mod_preds[6]), Metrics::rmse)
}

mod_evals12=function(split, model) {
  
  assessment_set= assessment(split) #get validation data
  
  holdout= assessment_set[,"abundance"] #select relevant column in dataframe
  
  mod_preds= predict(model, n.ahead=12)$pred
  
  rmse=pmap(list(holdout[12], mod_preds[12]), Metrics::rmse)
}

#for plotting frecasts
get_dat=function(split, model) {
  
  assessment_set= assessment(split) #get validation data
  
  holdout= assessment_set[,"abundance"]
  moon= assessment_set[,"newmoonnumber"]
  mod_preds=as.integer(predict(model, n.ahead=12, newxreg=assessment_set[,5:7])$pred)
  hp=cbind(moon, holdout, mod_preds)
}

###for plotting forecast evals

#RMSE~newmoon
get_dat1=function(split, evals_same, id) {
  
  assessment_set= assessment(split) #get validation data
  
  newmoon= assessment_set[,"newmoonnumber"][1]
  
  pred_score=evals_same
  id=id
  h=1
  hp=cbind(newmoon,id, pred_score, h)
}

get_dat6=function(split, evals_same6, id) {
  
  assessment_set= assessment(split) #get validation data
  
  newmoon= assessment_set[,"newmoonnumber"][6]
   pred_score=evals_same6
  id=id
  h=6
  hp=cbind(newmoon,id, pred_score, h)
}

get_dat12=function(split, evals_same12, id) {
  
  assessment_set= assessment(split) #get validation data
  
  newmoon= assessment_set[,"newmoonnumber"][12]
  pred_score=evals_same12
  id=id
  h=12
  hp=cbind(newmoon,id, pred_score, h)
}

