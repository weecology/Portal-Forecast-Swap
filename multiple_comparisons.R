#MULTIPLE COMPARISONS TEST####

#get p-value####

#calculate pvalue following: https://www.bmj.com/content/343/bmj.d2304 

get_pvalue=function(model) {

  coef_est=tscount::se(model)$est
  coef_ci_u=tscount::se(model)$ci[,"upper"]
  coef_ci_l=tscount::se(model)$ci[,"lower"]
  
  coef_se= (coef_ci_u - coef_ci_l)/ (2*1.96)
  
  coef_z=coef_est/coef_se
  
  coef_P = list(round(exp(-0.717 * coef_z - 0.416* ((coef_z)^2)), digits=4))
  
}

#add column of pvalues

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

#pairwise t-test####
pairwise.t.test(pbint_pval$raw_pvalue, pbint_pval$treatment, p.adjust.method="fdr") # <0.05
pairwise.t.test(pbb1_pval$raw_pvalue, pbb1_pval$treatment, p.adjust.method="fdr") #0.32
pairwise.t.test(pbb12_pval$raw_pvalue, pbb12_pval$treatment, p.adjust.method="fdr") #0.009
pairwise.t.test(pbtemp_pval$raw_pvalue, pbtemp_pval$treatment, p.adjust.method="fdr") #0.004
pairwise.t.test(pbwarmprec_pval$raw_pvalue, pbwarmprec_pval$treatment, p.adjust.method="fdr") #0.32
pairwise.t.test(pbcoolprec_pval$raw_pvalue, pbcoolprec_pval$treatment, p.adjust.method="fdr")  #<0.05

#add adjusted p-values (FDR)####
pbint_pval$FDR=p.adjust(pbint_pval$raw_pvalue, method="fdr")
pbb1_pval$FDR=p.adjust(pbb1_pval$raw_pvalue, method="fdr")
pbb12_pval$FDR=p.adjust(pbb12_pval$raw_pvalue, method="fdr")
pbtemp_pval$FDR=p.adjust(pbtemp_pval$raw_pvalue, method="fdr")
pbwarmprec_pval$FDR=p.adjust(pbwarmprec_pval$raw_pvalue, method="fdr")
pbcoolprec_pval$FDR=p.adjust(pbcoolprec_pval$raw_pvalue, method="fdr")
