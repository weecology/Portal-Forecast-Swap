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

control_b13_coefs=PPcontrol_dat$model%>%map(coef)%>%map_dbl(3)
exclosure_b13_coefs=PPexclosure_dat$model%>%map(coef)%>%map_dbl(3)

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

b13=cbind(control_b13_coefs, exclosure_b13_coefs)%>%as.data.frame%>%
  rename("control"="control_b13_coefs", "exclosure"="exclosure_b13_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "beta13")
coef_df=as.data.frame(list(ints, b1, b13, temps, warmprec, coolprec))%>%select(treatment, intercept, beta1, beta13, temp,cool_precip, warm_precip)

p1=ggviolin(coef_df, x="treatment", y="intercept", fill="treatment",
            palette=c("#69b3a2", "grey"),add="dotplot", add.params=list(fill="black"))+geom_hline(yintercept=0, lty=2)
p2=ggviolin(coef_df, x="treatment", y="beta1", fill="treatment",
            palette=c("#69b3a2", "grey"), add="dotplot", add.params=list(fill="black"))+geom_hline(yintercept=0, lty=2)
p3=ggviolin(coef_df, x="treatment", y="beta13", fill="treatment",
            palette=c("#69b3a2", "grey"),  add="dotplot", add.params=list(fill="black"))+geom_hline(yintercept=0, lty=2)
p4=ggviolin(coef_df, x="treatment", y="temp", fill="treatment",
            palette=c("#69b3a2", "grey"),  add="dotplot", add.params=list(fill="black"))+geom_hline(yintercept=0, lty=2)+ylab("mean temperature (lag=1)")
p5=ggviolin(coef_df, x="treatment", y="warm_precip", fill="treatment",
            palette=c("#69b3a2", "grey"),  add="dotplot", add.params=list(fill="black"))+ylab("warm precipitation")+geom_hline(yintercept=0, lty=2)
p6=ggviolin(coef_df, x="treatment", y="cool_precip", fill="treatment",
            palette=c("#69b3a2", "grey"),add="dotplot", add.params=list(fill="black"))+ylab("cool precipitation")+geom_hline(yintercept=0, lty=2)

coefs_plot=ggarrange(p1,p2,p3,p4,p5,p6, common.legend = T)
annotate_figure(coefs_plot, top=text_grob("C. penicillatus", face="bold", size=14))

#coefficient overlap####
int_ov=ggdensity(coef_df, x="intercept", color="treatment",fill="treatment", 
                 palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                 main="intercept")+geom_vline(xintercept=0)

ar1_ov=ggdensity(coef_df, x="beta1", color="treatment",fill="treatment", 
                 palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                 main="AR (1)")+geom_vline(xintercept=0)

ar13_ov=ggdensity(coef_df, x="beta13", color="treatment",fill="treatment", 
                  palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                  main="AR(13)")+geom_vline(xintercept=0)

temp_ov=ggdensity(coef_df, x="temp", color="treatment",fill="treatment", 
                  palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                  main="mean temperature (lag=1)")+geom_vline(xintercept=0)
wprec_ov=ggdensity(coef_df, x="warm_precip", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                   main="warm precipitation")+geom_vline(xintercept=0)
cprec_ov=ggdensity(coef_df, x="cool_precip", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                   main="cool precipitation")+geom_vline(xintercept=0)

coef_ov=ggarrange(int_ov, ar1_ov, ar13_ov, temp_ov, wprec_ov, cprec_ov, common.legend = T)
annotate_figure(coef_ov, top=text_grob("C. penicillatus", face="bold", size=14))

##forecasts####

#h=12
#control-control
cont_preds_same=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$preds_same), get_dat_same)

#control dat-exclosure mod
cont_preds_switch=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$preds_switch), get_dat_switch)

m=rep(seq(1:47), each=13) #no.of splits
code1="same"
code2="switched"

preds_cont_same=do.call(rbind.data.frame, cont_preds_same)
preds_cont_same=cbind(preds_cont_same, m, code1)
preds_cont_switch=do.call(rbind.data.frame, cont_preds_switch)
preds_cont_switch=cbind(preds_cont_switch, m, code2)

preds_plot1=ggplot()+
  geom_line(data=ppcont_covs_ro, aes(y=abundance, x=newmoonnumber), size=0.75)+
#  geom_point(data=preds_cont_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4, pch=19)+
 # geom_point(data=preds_cont_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4, pch=19)+
  geom_line(data=preds_cont_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4)+
  geom_line(data=preds_cont_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4)+
  theme_classic()+labs(colour="configuration")+
  ggtitle("PP control")+
  scale_colour_manual(values=c(same="darkblue", switched="darkred"), labels=c("matched", "mismatched"))

preds_plot1

#exclosure
excl_preds_same=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$preds_same), get_dat_same)

#exclosure dat-controlmod
excl_preds_switch=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$preds_switch), get_dat_switch)

preds_excl_same=do.call(rbind.data.frame, excl_preds_same)
preds_excl_same=cbind(preds_excl_same, m, code1)
preds_excl_switch=do.call(rbind.data.frame, excl_preds_switch)
preds_excl_switch=cbind(preds_excl_switch, m, code2)

preds_plot2=ggplot()+geom_line(data=ppexcl_covs_ro, aes(y=abundance, x=newmoonnumber), size=0.75)+
 # geom_point(data=preds_excl_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4, pch=19)+
#  geom_point(data=preds_excl_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4, pch=19)+
  geom_line(data=preds_excl_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4)+
  geom_line(data=preds_excl_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4)+
  theme_classic()+labs(colour="configuration")+
  ggtitle("PP exclosure")+
  scale_colour_manual(values=c(same="darkblue", switched="darkred"), labels=c("matched", "mismatched"))

preds_plot2

ggarrange(preds_plot1, preds_plot2, common.legend = T, ncol=1, nrow=2)

##forecast evals####

###RMSE differences####

cont_evals1_diff=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_same, PPcontrol_dat$evals_switch ,PPcontrol_dat$id), get_evals1_diff)
evals1_cont_diff=do.call(rbind.data.frame, cont_evals1_diff)%>%mutate(plot="control")

excl_evals1_diff=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_same, PPexclosure_dat$evals_switch ,PPexclosure_dat$id), get_evals1_diff)
evals1_excl_diff=do.call(rbind.data.frame, excl_evals1_diff)%>%mutate(plot="exclosure")

evals1=rbind(evals1_cont_diff, evals1_excl_diff)
evals1$newmoon=as.integer(evals1$newmoon)
evals1$h=as.integer(evals1$h)
evals1$score_same=as.numeric(evals1$score_same)
evals1$score_switch=as.numeric(evals1$score_switch)
evals1$score_diff=as.numeric(evals1$score_diff)

str(evals1)

#h=6
cont_evals6_diff=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_same, PPcontrol_dat$evals_switch ,PPcontrol_dat$id), get_evals6_diff)
evals6_cont_diff=do.call(rbind.data.frame, cont_evals6_diff)%>%mutate(plot="control")

excl_evals6_diff=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_same, PPexclosure_dat$evals_switch ,PPexclosure_dat$id), get_evals6_diff)
evals6_excl_diff=do.call(rbind.data.frame, excl_evals6_diff)%>%mutate(plot="exclosure")

evals6=rbind(evals6_cont_diff, evals6_excl_diff)
evals6$newmoon=as.integer(evals6$newmoon)
evals6$h=as.integer(evals6$h)
evals6$score_same=as.numeric(evals6$score_same)
evals6$score_switch=as.numeric(evals6$score_switch)
evals6$score_diff=as.numeric(evals6$score_diff)

str(evals6)

#h=12
cont_evals13_diff=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_same, PPcontrol_dat$evals_switch ,PPcontrol_dat$id), get_evals13_diff)
evals13_cont_diff=do.call(rbind.data.frame, cont_evals13_diff)%>%mutate(plot="control")

excl_evals13_diff=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_same, PPexclosure_dat$evals_switch ,PPexclosure_dat$id), get_evals13_diff)
evals13_excl_diff=do.call(rbind.data.frame, excl_evals13_diff)%>%mutate(plot="exclosure")

evals13=rbind(evals13_cont_diff, evals13_excl_diff)
evals13$newmoon=as.integer(evals13$newmoon)
evals13$h=as.integer(evals13$h)
evals13$score_same=as.numeric(evals13$score_same)
evals13$score_switch=as.numeric(evals13$score_switch)
evals13$score_diff=as.numeric(evals13$score_diff)

##plot###

pp_h1=ggdensity(evals1, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-25, 25)+ ylim(-0.01, 0.5)+
  annotate("segment", x=-5, y= 0.45, xend=-23, yend=0.45, col="black", size=1, arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=5, y= 0.45, xend=23, yend=0.45, col="black", size=1, arrow=arrow(length=unit(0.3, "cm")))+
  annotate("text", x=-15, y=0.40, label="matched better fit")+
  annotate("text", x=13, y=0.40, label="mismatched better fit")

pp_h6=ggdensity(evals6, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), rug=T, add="mean", xlab=F,size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-25, 25)+ylim(-0.01, 0.15)

pp_h13=ggdensity(evals13, x="score_diff", color="plot",fill="plot", size=1,
                 palette=c("#69b3a2", "grey"), rug=T, add="mean", xlab="RMSE difference (matched-mismatched)",
                 main="h=13")+geom_vline(xintercept=0)+xlim(-25, 25)+ ylim(-0.01, 0.15)

pph=ggarrange(pp_h1,pp_h6,pp_h13, common.legend = T, nrow=3)
annotate_figure(pph, top = text_grob("C. penicillatus", face = "bold", size = 14))


###Brier scores####

#control###

#combine predictions on same and switched models for control data
pp_preds_control=left_join(preds_cont_same, preds_cont_switch, by=c("moon", "holdout", "m"))

#calculate brier score for same models
pp_brier_cont1=scoring(pred=pp_preds_control$preds_same, response=pp_preds_control$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine same model predictions and brier score dataframes
pp_brier_cont_same=cbind(pp_preds_control, pp_brier_cont1)%>%rename(quadratic_same=quadratic)

#calculate brier score for switched models
pp_brier_cont2=scoring(pred=pp_preds_control$preds_switch, response=pp_preds_control$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine switched model predictions and brier score dataframes
pp_brier_cont_switch=cbind(pp_preds_control, pp_brier_cont2)%>%rename(quadratic_switch=quadratic)

#calculate brier score differences
pp_brier_control=left_join(pp_brier_cont_same, pp_brier_cont_switch)%>%mutate(treatment="control",brier_diff=quadratic_same-quadratic_switch)

#select relevant rows with horizons 1,6 13
pp_brier_control1=pp_brier_control%>%group_by(m)%>%slice(which(row_number()%in%c(1,6,13)))

#subdivide into each horizon for easier plotting
ppb1=pp_brier_control1%>%group_by(m)%>%slice_head(n=1)
ppb6=pp_brier_control1%>%group_by(m)%>%slice(which(row_number()==2))%>%mutate(horizon="6")
ppb13=pp_brier_control1%>%group_by(m)%>%slice(which(row_number()==3))%>%mutate(horizon="13")

#combine all 3 horizons
ppb=rbind(ppb1,ppb6,ppb13)

#exclosure###

#combine predictions on same and switched models for exclosure data
pp_preds_exclosure=left_join(preds_excl_same, preds_excl_switch, by=c("moon", "holdout", "m"))

#calculate brier score for same models
pp_brier_excl1=scoring(pred=pp_preds_exclosure$preds_same, response=pp_preds_exclosure$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine same model predictions and brier score dataframes
pp_brier_excl_same=cbind(pp_preds_exclosure, pp_brier_excl1)%>%rename(quadratic_same=quadratic)

#calculate brier score for switched models
pp_brier_excl2=scoring(pred=pp_preds_exclosure$preds_switch, response=pp_preds_exclosure$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine switched model predictions and brier score dataframes
pp_brier_excl_switch=cbind(pp_preds_exclosure, pp_brier_excl2)%>%rename(quadratic_switch=quadratic)

#calculate brier score differences
pp_brier_exclosure=left_join(pp_brier_excl_same, pp_brier_excl_switch)%>%mutate(treatment="exclosure",brier_diff=quadratic_same-quadratic_switch)

#select relevant rows with horizons 1,6 13
pp_brier_exclosure1=pp_brier_exclosure%>%group_by(m)%>%slice(which(row_number()%in%c(1,6,13)))

#subdivide into each horizon for easier plotting
ppbx1=pp_brier_exclosure1%>%group_by(m)%>%slice_head(n=1)%>%mutate(horizon="1")
ppbx6=pp_brier_exclosure1%>%group_by(m)%>%slice(which(row_number()==2))%>%mutate(horizon="6")
ppbx13=pp_brier_exclosure1%>%group_by(m)%>%slice(which(row_number()==3))%>%mutate(horizon="13")

#combine all 3 horizons
ppbx=rbind(ppbx1,ppbx6, ppbx13)

#combine control and exclosure and filter out per horizon for plotting
pp_briers1=rbind(ppb,ppbx)%>%filter(horizon==1)
pp_briers6=rbind(ppb,ppbx)%>%filter(horizon==6)
pp_briers13=rbind(ppb,ppbx)%>%filter(horizon==13)

pbp1=ggdensity(pp_briers1, x="brier_diff", color="treatment",fill="treatment", 
               palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
               main="h=1")+geom_vline(xintercept=0)+xlim(-0.30, .30)

pbp2=ggdensity(pp_briers6, x="brier_diff", color="treatment",fill="treatment", 
               palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
               main="h=6")+geom_vline(xintercept=0)+xlim(-0.30, .30)

pbp3=ggdensity(pp_briers13, x="brier_diff", color="treatment",fill="treatment", 
               palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
               main="h=13")+geom_vline(xintercept=0)+xlim(-0.30, .30)

ppbh=ggarrange(pbp1,pbp2,pbp3, common.legend = T, nrow=3)
annotate_figure(ppbh, top = text_grob("C. penicillatus", face = "bold", size = 14))

#####extra plot; unused###
###RMSE~moon###

#control

#h=1#
evals1_control=evals1%>%filter(plot=="control")
evals6_control=evals6%>%filter(plot=="control")
evals13_control=evals13%>%filter(plot=="control")

evals1_exclosure=evals1%>%filter(plot=="exclosure")
evals6_exclosure=evals6%>%filter(plot=="exclosure")
evals13_exclosure=evals13%>%filter(plot=="exclosure")

evals_moon1_cont=ggplot()+
  geom_line(data=evals1_control, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=evals1_control, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PP control (h=1)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

evals_moon1_cont

evals_moon6_cont=ggplot()+
  geom_line(data=evals6_control, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=evals6_control, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PP control (h=6)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

evals_moon6_cont

evals_moon13_cont=ggplot()+
  geom_line(data=evals13_control, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=evals13_control, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PP control (h=13)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

evals_moon13_cont

ggarrange(evals_moon1_cont, evals_moon6_cont, evals_moon13_cont, common.legend = T, nrow=1, ncol=3)

#exclosure
evals_moon1_excl=ggplot()+
  geom_line(data=evals1_exclosure, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=evals1_exclosure, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PP exclosure (h=1)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

evals_moon1_excl

evals_moon6_excl=ggplot()+
  geom_line(data=evals6_exclosure, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=evals6_exclosure, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PP exclosure (h=6)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

evals_moon6_excl

evals_moon13_excl=ggplot()+
  geom_line(data=evals12_exclosure, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=evals12_exclosure, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PP exclosure (h=12)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

evals_moon12_excl

ggarrange(evals_moon1_excl, evals_moon6_excl, evals_moon12_excl, common.legend = T, nrow=1, ncol=3)