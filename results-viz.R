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
            palette=c("#69b3a2", "grey"), add="jitter", add.params=list(fill="white"))
p2=ggviolin(coef_df, x="treatment", y="beta1", fill="treatment",
            palette=c("#69b3a2", "grey"), add="jitter", add.params=list(fill="white"))
p3=ggviolin(coef_df, x="treatment", y="beta13", fill="treatment",
            palette=c("#69b3a2", "grey"), add="jitter", add.params=list(fill="white"))
p4=ggviolin(coef_df, x="treatment", y="temp", fill="treatment",
            palette=c("#69b3a2", "grey"), add="jitter", add.params=list(fill="white"))+ylab("mean temperature (lag=1)")
p5=ggviolin(coef_df, x="treatment", y="warm_precip", fill="treatment",
            palette=c("#69b3a2", "grey"), add="jitter", add.params=list(fill="white"))+ylab("warm precipitation")
p6=ggviolin(coef_df, x="treatment", y="cool_precip", fill="treatment",
            palette=c("#69b3a2", "grey"), add="jitter", add.params=list(fill="white"))+ylab("cool precipitation")

coefs_plot=ggarrange(p1,p2,p3,p4,p5,p6, common.legend = T)
annotate_figure(coefs_plot, top=text_grob("PP coefficients", face="bold", size=14))

##forecasts####

#h=12
#control-control
cont_preds_same=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$preds_same), get_dat_same)

#control dat-exclosure mod
cont_preds_switch=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$preds_switch), get_dat_switch)

m=rep(seq(1:39), each=12) #no.of splits
code1="same"
code2="switched"

preds_cont_same=do.call(rbind.data.frame, cont_preds_same)
preds_cont_same=cbind(preds_cont_same, m, code1)
preds_cont_switch=do.call(rbind.data.frame, cont_preds_switch)
preds_cont_switch=cbind(preds_cont_switch, m, code2)

preds_plot1=ggplot()+
  geom_line(data=pp_datc, aes(y=abundance, x=newmoonnumber))+
  geom_point(data=preds_cont_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4, pch=19)+
  geom_point(data=preds_cont_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4, pch=19)+
  geom_line(data=preds_cont_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4)+
  geom_line(data=preds_cont_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4)+
  theme_classic()+labs(colour="configuration")+
  ggtitle("PP control")+
  scale_colour_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

preds_plot1

#exclosure
excl_preds_same=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$preds_same), get_dat_same)

#exclosure dat-controlmod
excl_preds_switch=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$preds_switch), get_dat_switch)

preds_excl_same=do.call(rbind.data.frame, excl_preds_same)
preds_excl_same=cbind(preds_excl_same, m, code1)
preds_excl_switch=do.call(rbind.data.frame, excl_preds_switch)
preds_excl_switch=cbind(preds_excl_switch, m, code2)

preds_plot2=ggplot()+geom_line(data=pp_date, aes(y=abundance, x=newmoonnumber))+
  geom_point(data=preds_excl_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4, pch=19)+
  geom_point(data=preds_excl_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4, pch=19)+
  geom_line(data=preds_excl_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4)+
  geom_line(data=preds_excl_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4)+
  theme_classic()+labs(colour="configuration")+
  ggtitle("PP exclosure")+
  scale_colour_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

preds_plot2

ggarrange(preds_plot1, preds_plot2, common.legend = T)

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

#h=1
pp_h1=ggplot(evals1, aes(score_diff, colour = plot, fill=plot))+
  geom_histogram(alpha=0.4, position="identity")+theme_classic()+xlab("RMSE difference (same-switched)")+
  scale_fill_manual(values=c(control="#69b3a2", exclosure="grey"))+
  scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))+
  geom_vline(xintercept=0, lty=2)+ggtitle("PP (h=1)")+facet_wrap(~plot)

pp_h1

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

pp_h6=ggplot(evals6, aes(score_diff, colour = plot, fill=plot))+
  geom_histogram(alpha=0.4, position="identity")+theme_classic()+xlab("RMSE difference (same-switched)")+
  scale_fill_manual(values=c(control="#69b3a2", exclosure="grey"))+
  scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))+
  geom_vline(xintercept=0, lty=2)+ggtitle("PP (h=6)")+facet_wrap(~plot)

pp_h6

#h=12
cont_evals12_diff=pmap(list(PPcontrol_dat$splits,PPcontrol_dat$evals_same, PPcontrol_dat$evals_switch ,PPcontrol_dat$id), get_evals12_diff)
evals12_cont_diff=do.call(rbind.data.frame, cont_evals12_diff)%>%mutate(plot="control")

excl_evals12_diff=pmap(list(PPexclosure_dat$splits,PPexclosure_dat$evals_same, PPexclosure_dat$evals_switch ,PPexclosure_dat$id), get_evals12_diff)
evals12_excl_diff=do.call(rbind.data.frame, excl_evals12_diff)%>%mutate(plot="exclosure")

evals12=rbind(evals12_cont_diff, evals12_excl_diff)
evals12$newmoon=as.integer(evals12$newmoon)
evals12$h=as.integer(evals12$h)
evals12$score_same=as.numeric(evals12$score_same)
evals12$score_switch=as.numeric(evals12$score_switch)
evals12$score_diff=as.numeric(evals12$score_diff)

str(evals12)

pp_h12=ggplot(evals12, aes(score_diff, colour = plot, fill=plot))+
  geom_histogram(alpha=0.4, position="identity")+theme_classic()+xlab("RMSE difference (same-switched)")+
  scale_fill_manual(values=c(control="#69b3a2", exclosure="grey"))+
  scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))+
  geom_vline(xintercept=0, lty=2)+ggtitle("PP (h=12)")+facet_wrap(~plot)

pp_h12

ggarrange(pp_h1,pp_h6,pp_h12, common.legend = T)

###RMSE~moon#####

#control

#h=1#
evals1_control=evals1%>%filter(plot=="control")
evals6_control=evals6%>%filter(plot=="control")
evals12_control=evals12%>%filter(plot=="control")

evals1_exclosure=evals1%>%filter(plot=="exclosure")
evals6_exclosure=evals6%>%filter(plot=="exclosure")
evals12_exclosure=evals12%>%filter(plot=="exclosure")

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

evals_moon12_cont=ggplot()+
  geom_line(data=evals12_control, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=evals12_control, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PP control (h=12)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

evals_moon12_cont

ggarrange(evals_moon1_cont, evals_moon6_cont, evals_moon12_cont, common.legend = T, nrow=1, ncol=3)

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

evals_moon12_excl=ggplot()+
  geom_line(data=evals12_exclosure, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=evals12_exclosure, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PP exclosure (h=12)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

evals_moon12_excl

ggarrange(evals_moon1_excl, evals_moon6_excl, evals_moon12_excl, common.legend = T, nrow=1, ncol=3)