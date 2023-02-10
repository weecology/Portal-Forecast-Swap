#results viz####

require(dplyr)
require(tidyr)
require(ggplot2)
require(ggpubr)

##coefficients####
control_int_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(1)
exclosure_int_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(1)

control_b1_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(2)
exclosure_b1_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(2)

control_a12_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(3)
exclosure_a12_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(3)

control_temp1_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(4)
exclosure_temp1_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(4)

control_warmprec_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(5)
exclosure_warmprec_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(5)

control_coolprec_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(6)
exclosure_coolprec_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(6)

warmprec=cbind(control_warmprec_coefs, exclosure_warmprec_coefs)%>%as.data.frame%>%
  rename("control"="control_warmprec_coefs", "exclosure"="exclosure_warmprec_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "warm_precip")

coolprec=cbind(control_coolprec_coefs, exclosure_coolprec_coefs)%>%as.data.frame%>%
  rename("control"="control_coolprec_coefs", "exclosure"="exclosure_coolprec_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "cool_precip")

temps=cbind(control_temp1_coefs, exclosure_temp1_coefs)%>%as.data.frame%>%
  rename("control"="control_temp1_coefs", "exclosure"="exclosure_temp1_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "temp")

ints=cbind(control_int_coefs, exclosure_int_coefs)%>%as.data.frame%>%
  rename("control"="control_int_coefs", "exclosure"="exclosure_int_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "intercept")

b1=cbind(control_b1_coefs, exclosure_b1_coefs)%>%as.data.frame%>%
  rename("control"="control_b1_coefs", "exclosure"="exclosure_b1_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "beta1")

a12=cbind(control_a12_coefs, exclosure_a12_coefs)%>%as.data.frame%>%
  rename("control"="control_a12_coefs", "exclosure"="exclosure_a12_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "alpha12")

coef_df=as.data.frame(list(ints, b1, a12, temps, warmprec, coolprec))%>%select(treatment, intercept, beta1, alpha12, temp,cool_precip, warm_precip)

p1=ggplot(coef_df, aes(x=treatment, y=intercept, col=treatment))+geom_boxplot()+theme_classic()+ylab("mean coefficient estimate")+ggtitle("intercept")
p2=ggplot(coef_df, aes(x=treatment, y=beta1, col=treatment))+geom_boxplot()+theme_classic()+ylab("mean coefficient estimate")+ggtitle("beta1")
p3=ggplot(coef_df, aes(x=treatment, y=alpha12, col=treatment))+geom_boxplot()+theme_classic()+ylab("mean coefficient estimate")+ggtitle("alpha12")
p4=ggplot(coef_df, aes(x=treatment, y=temp, col=treatment))+geom_boxplot()+theme_classic()+ylab("mean coefficient estimate")+ggtitle("mean temperature (lag 1)")
p5=ggplot(coef_df, aes(x=treatment, y=warm_precip, col=treatment))+geom_boxplot()+theme_classic()+ylab("mean coefficient estimate")+ggtitle("warm precipitation")
p6=ggplot(coef_df, aes(x=treatment, y=cool_precip, col=treatment))+geom_boxplot()+theme_classic()+ylab("mean coefficient estimate")+ggtitle("cool precipitation")

coefs_plot=ggarrange(p1,p2,p3,p4,p5,p6, common.legend = T)
annotate_figure(coefs_plot, top=text_grob("PP coefficients", face="bold", size=14))

##forecasts####

#h=12
#control-control
cont_preds_same=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$model), get_dat)

#control dat-exclosure mod
cont_preds_switch=pmap(list(PPcontrol_dat$splits,PPexclosure_dat$model), get_dat)

m=rep(seq(1:39), each=12)
code1="same"
code2="switched"

preds_cont_same=do.call(rbind.data.frame, cont_preds_same)
preds_cont_same=cbind(preds_cont_same, m, code1)
preds_cont_switch=do.call(rbind.data.frame, cont_preds_switch)
preds_cont_switch=cbind(preds_cont_switch, m, code2)

preds_plot1=ggplot()+geom_line(data=pp_datc, aes(y=abundance, x=newmoonnumber))+
  geom_line(data=preds_cont_same, aes(y=mod_preds, x=moon, group=h,color="same"), alpha=0.4)+
  geom_line(data=preds_cont_switch, aes(y=mod_preds, x=moon, group=h, color="switched"), alpha=0.4)+
  theme_classic()+
  ggtitle("PP control")+
  scale_colour_manual(values=c(same="blue", switched="red"), labels=c("same", "switched"))

preds_plot1

#exclosure
excl_preds_same=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$model), get_dat)

#exclosure dat-controlmod
excl_preds_switch=pmap(list(PPexclosure_dat$splits,PPcontrol_dat$model), get_dat)

preds_excl_same=do.call(rbind.data.frame, excl_preds_same)
preds_excl_same=cbind(preds_excl_same, m, code1)
preds_excl_switch=do.call(rbind.data.frame, excl_preds_switch)
preds_excl_switch=cbind(preds_excl_switch, m, code2)

preds_plot2=ggplot()+geom_line(data=pp_date, aes(y=abundance, x=newmoonnumber))+
  geom_line(data=preds_excl_same, aes(y=mod_preds, x=moon, group=h,color="same"), alpha=0.4)+
  geom_line(data=preds_excl_switch, aes(y=mod_preds, x=moon, group=h, color="switched"), alpha=0.4)+
  theme_classic()+
  ggtitle("PP exclosure")+
  scale_colour_manual(values=c(same="blue", switched="red"), labels=c("same", "switched"))

preds_plot2

ggarrange(preds_plot1, preds_plot2, common.legend = T)

##forecast evals####

###RMSE differences####

#control
ppcontrol_evals1_same=as.vector(unlist(PPcontrol_dat$evals_same))
ppcontrol_evals1_switch=as.vector(unlist(PPcontrol_dat$evals_switch))

ppcontrol_evals6_same=as.vector(unlist(PPcontrol_dat$evals_same6))
ppcontrol_evals6_switch=as.vector(unlist(PPcontrol_dat$evals_switch6))

ppcontrol_evals12_same=as.vector(unlist(PPcontrol_dat$evals_same12))
ppcontrol_evals12_switch=as.vector(unlist(PPcontrol_dat$evals_switch12))

ppcontrol_evals1=as.data.frame(cbind(ppcontrol_evals1_same, ppcontrol_evals1_switch))
ppcontrol_evals6=as.data.frame(cbind(ppcontrol_evals6_same, ppcontrol_evals6_switch))
ppcontrol_evals12=as.data.frame(cbind(ppcontrol_evals12_same, ppcontrol_evals12_switch))

#exclosure
ppexclosure_evals1_same=as.vector(unlist(PPexclosure_dat$evals_same))
ppexclosure_evals1_switch=as.vector(unlist(PPexclosure_dat$evals_switch))

ppexclosure_evals6_same=as.vector(unlist(PPexclosure_dat$evals_same6))
ppexclosure_evals6_switch=as.vector(unlist(PPexclosure_dat$evals_switch6))

ppexclosure_evals12_same=as.vector(unlist(PPexclosure_dat$evals_same12))
ppexclosure_evals12_switch=as.vector(unlist(PPexclosure_dat$evals_switch12))

ppexclosure_evals1=as.data.frame(cbind(ppexclosure_evals1_same, ppexclosure_evals1_switch))
ppexclosure_evals6=as.data.frame(cbind(ppexclosure_evals6_same, ppexclosure_evals6_switch))
ppexclosure_evals12=as.data.frame(cbind(ppexclosure_evals12_same, ppexclosure_evals12_switch))

pp_evals1=cbind(ppcontrol_evals1, ppexclosure_evals1)%>%mutate(h="1")%>%
  rename("cont_dat_mod"="ppcontrol_evals1_same", "cont_dat_excl_mod"="ppcontrol_evals1_switch", 
         "excl_dat_mod"="ppexclosure_evals1_same",  "excl_dat_cont_mod"="ppexclosure_evals1_switch")

pp_evals6=cbind(ppcontrol_evals6, ppexclosure_evals6)%>%mutate(h="6")%>%
  rename("cont_dat_mod"="ppcontrol_evals6_same", "cont_dat_excl_mod"="ppcontrol_evals6_switch", 
         "excl_dat_mod"="ppexclosure_evals6_same",  "excl_dat_cont_mod"="ppexclosure_evals6_switch")

pp_evals12=cbind(ppcontrol_evals12, ppexclosure_evals12)%>%mutate(h="12")%>%
  rename("cont_dat_mod"="ppcontrol_evals12_same", "cont_dat_excl_mod"="ppcontrol_evals12_switch", 
         "excl_dat_mod"="ppexclosure_evals12_same",  "excl_dat_cont_mod"="ppexclosure_evals12_switch")

pp_evals=rbind(pp_evals1, pp_evals6, pp_evals12)

#for RMSE diff plots

#RMSE difference per plot
#h=1
ppevals_a=pp_evals%>%filter(h==1)%>%select(1:2, h)%>%
  mutate(plot="control", config_diff=cont_dat_mod-cont_dat_excl_mod)%>%
  select(h, plot, config_diff)

ppevals_b=pp_evals%>%filter(h==1)%>%select(3:4, h)%>%
  mutate(plot="exclosure", config_diff=excl_dat_mod-excl_dat_cont_mod)%>%
  select(h, plot, config_diff)

ppevals_c=rbind(ppevals_a, ppevals_b)

a1=ggplot(ppevals_c, aes(config_diff, colour = plot, fill=plot))+
  geom_histogram(alpha=0.4, position="identity")+theme_classic()+xlab("RMSE difference (same-switched)")+
  geom_vline(xintercept=0, lty=2)+ggtitle("PP (h=1)")+facet_wrap(~plot)

#h=6
ppevals6_a=pp_evals%>%filter(h==6)%>%select(1:2, h)%>%
  mutate(plot="control", config_diff=cont_dat_mod-cont_dat_excl_mod)%>%
  select(h, plot, config_diff)

ppevals6_b=pp_evals%>%filter(h==6)%>%select(3:4, h)%>%
  mutate(plot="exclosure", config_diff=excl_dat_mod-excl_dat_cont_mod)%>%
  select(h, plot, config_diff)

ppevals6_c=rbind(ppevals6_a, ppevals6_b)

b1=ggplot(ppevals6_c, aes(config_diff, colour = plot, fill=plot))+
  geom_histogram(alpha=0.4, position="identity")+theme_classic()+xlab("RMSE difference (same-switched)")+
  geom_vline(xintercept=0, lty=2)+ggtitle("PP (h=6)")+facet_wrap(~plot)

#h=12
ppevals12_a=pp_evals%>%filter(h==12)%>%select(1:2, h)%>%
  mutate(plot="control", config_diff=cont_dat_mod-cont_dat_excl_mod)%>%
  select(h, plot, config_diff)

ppevals12_b=pp_evals%>%filter(h==12)%>%select(3:4, h)%>%
  mutate(plot="exclosure", config_diff=excl_dat_mod-excl_dat_cont_mod)%>%
  select(h, plot, config_diff)

ppevals12_c=rbind(ppevals12_a, ppevals12_b)

c1=ggplot(ppevals12_c, aes(config_diff, colour = plot, fill=plot))+
  geom_histogram(alpha=0.4, position="identity")+theme_classic()+xlab("RMSE difference (same-switched)")+
  geom_vline(xintercept=0, lty=2)+ggtitle("PP (h=12)")+facet_wrap(~plot)

ggarrange(a1,b1,c1, common.legend = T)

#for RMSE~period and RMSE~h plots
ppevals=pp_evals%>%pivot_longer(cols=1:4, names_to = "config", values_to="RMSE")%>%
  mutate(configuration=case_when(
    config=="cont dat-mod" ~"same",
    config=="cont dat-excl mod" ~"switched",
    config=="excl dat-mod" ~"same",
    config=="excl dat-cont mod" ~"switched"
  ),
  plot=case_when( config=="cont dat-mod" ~"control",
                  config=="cont dat-excl mod" ~"control",
                  config=="excl dat-mod" ~"exclosure",
                  config=="excl dat-cont mod" ~"exclosure"))

ppevals_h1=ppevals%>%filter(m==1)
ppevals_h6=ppevals%>%filter(m==6)
ppevals_h12=ppevals%>%filter(m==12)

###RMSE~period#####

#h=1###
#control
pc1=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_same, PPcontrol_dat$id), get_dat1)
pc2=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_switch, PPcontrol_dat$id), get_dat1)

pc1_lst1=pc1%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="same")

pc1_lst2=pc2%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="switched")

pc1_lst1$newmoon=as.integer(pc1_lst1$newmoon)
pc1_lst1$pred_score=as.integer(pc1_lst1$pred_score)
pc1_lst1$h=as.integer(pc1_lst1$h)

pc1_lst2$newmoon=as.integer(pc1_lst2$newmoon)
pc1_lst2$pred_score=as.integer(pc1_lst2$pred_score)
pc1_lst2$h=as.integer(pc1_lst2$h)

lst1_c=rbind(pc1_lst1, pc1_lst2)
lst1_c=as.data.frame(lst1_c)%>%mutate(plot="control")
lst1_c$h=as.factor(lst1_c$h)

pc1_lst1$id=as.vector(unlist(pc1_lst1$id))
pc1_lst2$id=as.vector(unlist(pc1_lst2$id))

cont_period1=ggplot(lst1_c, aes(x=newmoon, y=pred_score, color=configuration))+
  geom_line(alpha=0.4)+theme_classic()+ylab("RMSE")+ggtitle("PP control (h=1)")

#exclosure
pe1=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_same, PPexclosure_dat$id), get_dat1)
pe2=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_switch, PPexclosure_dat$id), get_dat1)

pe1_lst1=pe1%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="same")

pe1_lst2=pe2%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="switched")

pe1_lst1$newmoon=as.integer(pe1_lst1$newmoon)
pe1_lst1$pred_score=as.integer(pe1_lst1$pred_score)
pe1_lst1$h=as.integer(pe1_lst1$h)

pe1_lst2$newmoon=as.integer(pe1_lst2$newmoon)
pe1_lst2$pred_score=as.integer(pe1_lst2$pred_score)
pe1_lst2$h=as.integer(pe1_lst2$h)

lst1_e=rbind(pe1_lst1, pe1_lst2)
lst1_e=as.data.frame(lst1_e)%>%mutate(plot="exclosure")
lst1_e$h=as.factor(lst1_e$h)

pe1_lst1$id=as.vector(unlist(pe1_lst1$id))
pe1_lst2$id=as.vector(unlist(pe1_lst2$id))

excl_period1=ggplot(lst1_e, aes(x=newmoon, y=pred_score, color=configuration))+
  geom_line(alpha=0.4)+theme_classic()+ylab("RMSE")+ggtitle("PP exclosure (h=1)")

h1_evals=rbind(lst1_c, lst1_e)

ggarrange(cont_period1, excl_period1, common.legend = T)

#h=6###
#control
pc6=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_same6, PPcontrol_dat$id), get_dat6)
pc7=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_switch6, PPcontrol_dat$id), get_dat6)

pc6_lst1=pc6%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="same")

pc6_lst2=pc7%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="switched")

pc6_lst1$newmoon=as.integer(pc6_lst1$newmoon)
pc6_lst1$pred_score=as.integer(pc6_lst1$pred_score)
pc6_lst1$h=as.integer(pc6_lst1$h)

pc6_lst2$newmoon=as.integer(pc6_lst2$newmoon)
pc6_lst2$pred_score=as.integer(pc6_lst2$pred_score)
pc6_lst2$h=as.integer(pc6_lst2$h)

pc6_lst1$id=as.vector(unlist(pc6_lst1$id))
pc6_lst2$id=as.vector(unlist(pc6_lst2$id))

lst6_c=rbind(pc6_lst1, pc6_lst2)%>%mutate(plot="control")
lst6_c$h=as.factor(lst6_c$h)

cont_period6=ggplot(lst6_c, aes(x=newmoon, y=pred_score, color=configuration))+
  geom_line(alpha=0.4)+theme_classic()+ylab("RMSE")+ggtitle("PP control (h=6)")

cont_period6

#exclosure
pe6=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_same6, PPexclosure_dat$id), get_dat6)
pe7=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_switch6, PPexclosure_dat$id), get_dat6)

pe6_lst1=pe6%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="same")

pe6_lst2=pe7%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="switched")

pe6_lst1$newmoon=as.integer(pe6_lst1$newmoon)
pe6_lst1$pred_score=as.integer(pe6_lst1$pred_score)
pe6_lst1$h=as.integer(pe6_lst1$h)

pe6_lst2$newmoon=as.integer(pe6_lst2$newmoon)
pe6_lst2$pred_score=as.integer(pe6_lst2$pred_score)
pe6_lst2$h=as.integer(pe6_lst2$h)

lst6_e=rbind(pe6_lst1, pe6_lst2)%>%mutate(plot="exclosure")
lst6_e$h=as.factor(lst6_e$h)
pe6_lst1$id=as.vector(unlist(pe6_lst1$id))
pe6_lst2$id=as.vector(unlist(pe6_lst2$id))

exc_period6=ggplot(lst6_e, aes(x=newmoon, y=pred_score, color=configuration))+
  geom_line(alpha=0.4)+theme_classic()+ylab("RMSE")+ggtitle("PP exclosure (h=6)")

exc_period6

ggarrange(cont_period6,exc_period6, common.legend = T)

#h=12###
#control 
pc12=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_same12, PPcontrol_dat$id), get_dat12)
pc13=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_switch12, PPcontrol_dat$id), get_dat12)

pc12_lst1=pc12%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="same")

pc12_lst2=pc13%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="switched")

pc12_lst1$newmoon=as.integer(pc12_lst1$newmoon)
pc12_lst1$pred_score=as.integer(pc12_lst1$pred_score)
pc12_lst1$h=as.integer(pc12_lst1$h)

pc12_lst2$newmoon=as.integer(pc12_lst2$newmoon)
pc12_lst2$pred_score=as.integer(pc12_lst2$pred_score)
pc12_lst2$h=as.integer(pc12_lst2$h)

pc12_lst1$id=as.vector(unlist(pc12_lst1$id))
pc12_lst2$id=as.vector(unlist(pc12_lst2$id))

lst12_c=rbind(pc12_lst1, pc12_lst2)%>%mutate(plot="control")
lst12_c$h=as.factor(lst12_c$h)

cont_period12=ggplot(lst12_c, aes(x=newmoon, y=pred_score, color=configuration))+
  geom_line(alpha=0.4)+theme_classic()+ylab("RMSE")+ggtitle("PP control(h=12)")

cont_period12


#exclosure
pe12=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_same12, PPexclosure_dat$id), get_dat12)
pe13=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_switch12, PPexclosure_dat$id), get_dat12)

pe12_lst1=pe12%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="same")

pe12_lst2=pe13%>%map_dfr(~as.data.frame(as.matrix(.)))%>%
  rename_with(~sub("V", "column",.))%>%mutate(configuration="switched")

pe12_lst1$newmoon=as.integer(pe12_lst1$newmoon)
pe12_lst1$pred_score=as.integer(pe12_lst1$pred_score)
pe12_lst1$h=as.integer(pe12_lst1$h)

pe12_lst2$newmoon=as.integer(pe12_lst2$newmoon)
pe12_lst2$pred_score=as.integer(pe12_lst2$pred_score)
pe12_lst2$h=as.integer(pe12_lst2$h)

lst12_e=rbind(pe12_lst1, pe12_lst2)%>%mutate(plot="exclosure")
lst12_e$h=as.factor(lst12_e$h)
pe12_lst1$id=as.vector(unlist(pe12_lst1$id))
pe12_lst2$id=as.vector(unlist(pe12_lst2$id))

exc_period12=ggplot(lst12_e, aes(x=newmoon, y=pred_score, color=configuration))+
  geom_line(alpha=0.4)+theme_classic()+ylab("RMSE")+ggtitle("PP exclosure (h=12)")

exc_period12

ggarrange(cont_period12,exc_period12, common.legend = T)
