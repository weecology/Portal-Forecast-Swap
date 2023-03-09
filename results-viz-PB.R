#results viz for PBs (data:2000-2010)####

require(dplyr)
require(tidyr)
require(ggplot2)
require(ggpubr)

##coefficients####
pbcontrol_int_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(1)
pbexclosure_int_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(1)

pbcontrol_b1_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(2)
pbexclosure_b1_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(2)

pbcontrol_b13_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(3)
pbexclosure_b13_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(3)

pbcontrol_temp1_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(4)
pbexclosure_temp1_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(4)

pbcontrol_warmprec_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(5)
pbexclosure_warmprec_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(5)

pbcontrol_coolprec_coefs=PBcontrol_dat$model%>%map(coef)%>%map_dbl(6)
pbexclosure_coolprec_coefs=PBexclosure_dat$model%>%map(coef)%>%map_dbl(6)

pbwarmprec=cbind(pbcontrol_warmprec_coefs, pbexclosure_warmprec_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_warmprec_coefs", "exclosure"="pbexclosure_warmprec_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "warm_precip")

pbcoolprec=cbind(pbcontrol_coolprec_coefs, pbexclosure_coolprec_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_coolprec_coefs", "exclosure"="pbexclosure_coolprec_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "cool_precip")

pbtemps=cbind(pbcontrol_temp1_coefs, pbexclosure_temp1_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_temp1_coefs", "exclosure"="pbexclosure_temp1_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "temp")

pbints=cbind(pbcontrol_int_coefs, pbexclosure_int_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_int_coefs", "exclosure"="pbexclosure_int_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "intercept")

pbb1=cbind(pbcontrol_b1_coefs, pbexclosure_b1_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_b1_coefs", "exclosure"="pbexclosure_b1_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "beta1")

pbb13=cbind(pbcontrol_b13_coefs, pbexclosure_b13_coefs)%>%as.data.frame%>%
  rename("control"="pbcontrol_b13_coefs", "exclosure"="pbexclosure_b13_coefs")%>%
  pivot_longer(cols=c(1:2),names_to="treatment", values_to = "beta13")

coef_df_PB=as.data.frame(list(pbints, pbb1, pbb13, pbtemps, pbwarmprec, pbcoolprec))%>%
  select(treatment, intercept, beta1, beta13, temp,cool_precip, warm_precip)

pb1=ggviolin(coef_df_PB, x="treatment", y="intercept", fill="treatment",
             palette=c("#69b3a2", "grey"),add="dotplot", add.params=list(fill="black"))+geom_hline(yintercept=0, lty=2)
pb2=ggviolin(coef_df_PB, x="treatment", y="beta1", fill="treatment",
             palette=c("#69b3a2", "grey"),add="dotplot", add.params=list(fill="black"))+geom_hline(yintercept=0, lty=2)
pb3=ggviolin(coef_df_PB, x="treatment", y="beta13", fill="treatment",
             palette=c("#69b3a2", "grey"),add="dotplot", add.params=list(fill="black"))+geom_hline(yintercept=0, lty=2)
pb4=ggviolin(coef_df_PB, x="treatment", y="temp", fill="treatment",
             palette=c("#69b3a2", "grey"),add="dotplot", add.params=list(fill="black"))+geom_hline(yintercept=0, lty=2)+ylab("mean temperature (lag=1)")
pb5=ggviolin(coef_df_PB, x="treatment", y="warm_precip", fill="treatment",
             palette=c("#69b3a2", "grey"),add="dotplot", add.params=list(fill="black"))+geom_hline(yintercept=0, lty=2)+ylab("warm precipitation")
pb6=ggviolin(coef_df_PB, x="treatment", y="cool_precip", fill="treatment",
             palette=c("#69b3a2", "grey"),add="dotplot", add.params=list(fill="black"))+geom_hline(yintercept=0, lty=2)+ylab("cool precipitation")

coefs_plot_PB=ggarrange(pb1,pb2,pb3,pb4,pb5,pb6, common.legend = T)
annotate_figure(coefs_plot_PB, top=text_grob("C. baileyi", face="bold", size=14))

##forecasts####

#h=12
#control-control
pbcont_preds_same=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$preds_same), get_dat_same)

#control dat-exclosure mod
pbcont_preds_switch=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$preds_switch), get_dat_switch)

m=rep(seq(1:50), each=12)
code1="same"
code2="switched"

PBpreds_cont_same=do.call(rbind.data.frame, pbcont_preds_same)
PBpreds_cont_same=cbind(PBpreds_cont_same, m, code1)
PBpreds_cont_switch=do.call(rbind.data.frame, pbcont_preds_switch)
PBpreds_cont_switch=cbind(PBpreds_cont_switch, m, code2)

pbpreds_plot1=ggplot()+
  geom_line(data=pb_datc, aes(y=abundance, x=newmoonnumber), size=0.75)+
  #geom_point(data=PBpreds_cont_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4, pch=19)+
  #geom_point(data=PBpreds_cont_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4, pch=19)+
  geom_line(data=PBpreds_cont_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4)+
  geom_line(data=PBpreds_cont_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4)+
  theme_classic()+labs(colour="configuration")+
  ggtitle("PB control")+
  scale_colour_manual(values=c(same="darkblue", switched="darkred"), labels=c("matched", "mismatched"))

pbpreds_plot1

#exclosure
pbexcl_preds_same=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$preds_same), get_dat_same)

#exclosure dat-controlmod
pbexcl_preds_switch=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$preds_switch), get_dat_switch)

PBpreds_excl_same=do.call(rbind.data.frame, pbexcl_preds_same)
PBpreds_excl_same=cbind(PBpreds_excl_same, m, code1)
PBpreds_excl_switch=do.call(rbind.data.frame, pbexcl_preds_switch)
PBpreds_excl_switch=cbind(PBpreds_excl_switch, m, code2)

pbpreds_plot2=ggplot()+geom_line(data=pb_date, aes(y=abundance, x=newmoonnumber), size=0.75)+
 # geom_point(data=PBpreds_excl_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4, pch=19)+
  #geom_point(data=PBpreds_excl_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4, pch=19)+
  geom_line(data=PBpreds_excl_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4)+
  geom_line(data=PBpreds_excl_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4)+
  theme_classic()+labs(colour="configuration")+
  ggtitle("PB exclosure")+
  scale_colour_manual(values=c(same="darkblue", switched="darkred"), labels=c("matched", "mismatched"))

pbpreds_plot2

ggarrange(pbpreds_plot1, pbpreds_plot2, common.legend = T, nrow=2, ncol=1)

##forecast evals####

###RMSE differences####

pbcont_evals1_diff=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$evals_same, PBcontrol_dat$evals_switch ,PBcontrol_dat$id), get_evals1_diff)
pbevals1_cont_diff=do.call(rbind.data.frame, pbcont_evals1_diff)%>%mutate(plot="control")

pbexcl_evals1_diff=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$evals_same, PBexclosure_dat$evals_switch ,PBexclosure_dat$id), get_evals1_diff)
pbevals1_excl_diff=do.call(rbind.data.frame, pbexcl_evals1_diff)%>%mutate(plot="exclosure")

pbevals1=rbind(pbevals1_cont_diff, pbevals1_excl_diff)
pbevals1$newmoon=as.integer(pbevals1$newmoon)
pbevals1$h=as.integer(epbvals1$h)
pbevals1$score_same=as.numeric(pbevals1$score_same)
pbevals1$score_switch=as.numeric(pbevals1$score_switch)
pbevals1$score_diff=as.numeric(pbevals1$score_diff)

str(pbevals1)

#h=1
pb_h1=ggplot(pbevals1, aes(score_diff, colour = plot, fill=plot))+
  geom_histogram(alpha=0.4, position="identity")+theme_classic()+xlab("RMSE difference (same-switched)")+
  scale_fill_manual(values=c(control="#69b3a2", exclosure="grey"))+
  scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))+
  geom_vline(xintercept=0, lty=2)+ggtitle("PB (h=1)")+facet_wrap(~plot)

pb_h1

#h=6
pbcont_evals6_diff=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$evals_same, PBcontrol_dat$evals_switch ,PBcontrol_dat$id), get_evals6_diff)
pbevals6_cont_diff=do.call(rbind.data.frame, pbcont_evals6_diff)%>%mutate(plot="control")

pbexcl_evals6_diff=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$evals_same, PBexclosure_dat$evals_switch ,PBexclosure_dat$id), get_evals6_diff)
pbevals6_excl_diff=do.call(rbind.data.frame, pbexcl_evals6_diff)%>%mutate(plot="exclosure")

pbevals6=rbind(pbevals6_cont_diff, pbevals6_excl_diff)
pbevals6$newmoon=as.integer(pbevals6$newmoon)
pbevals6$h=as.integer(pbevals6$h)
pbevals6$score_same=as.numeric(pbevals6$score_same)
pbevals6$score_switch=as.numeric(pbevals6$score_switch)
pbevals6$score_diff=as.numeric(pbevals6$score_diff)

str(pbevals6)

pb_h6=ggplot(pbevals6, aes(score_diff, colour = plot, fill=plot))+
  geom_histogram(alpha=0.4, position="identity")+theme_classic()+xlab("RMSE difference (same-switched)")+
  scale_fill_manual(values=c(control="#69b3a2", exclosure="grey"))+
  scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))+
  geom_vline(xintercept=0, lty=2)+ggtitle("PB (h=6)")+facet_wrap(~plot)

pb_h6

#h=12
pbcont_evals12_diff=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$evals_same, PBcontrol_dat$evals_switch ,PBcontrol_dat$id), get_evals12_diff)
pbevals12_cont_diff=do.call(rbind.data.frame, pbcont_evals12_diff)%>%mutate(plot="control")

pbexcl_evals12_diff=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$evals_same, PBexclosure_dat$evals_switch ,PBexclosure_dat$id), get_evals12_diff)
pbevals12_excl_diff=do.call(rbind.data.frame, pbexcl_evals12_diff)%>%mutate(plot="exclosure")

pbevals12=rbind(pbevals12_cont_diff, pbevals12_excl_diff)
pbevals12$newmoon=as.integer(pbevals12$newmoon)
pbevals12$h=as.integer(pbevals12$h)
pbevals12$score_same=as.numeric(pbevals12$score_same)
pbevals12$score_switch=as.numeric(pbevals12$score_switch)
pbevals12$score_diff=as.numeric(pbevals12$score_diff)

str(pbevals12)

pb_h12=ggplot(pbevals12, aes(score_diff, colour = plot, fill=plot))+
  geom_histogram(alpha=0.4, position="identity")+theme_classic()+xlab("RMSE difference (same-switched)")+
  scale_fill_manual(values=c(control="#69b3a2", exclosure="grey"))+
  scale_color_manual(values=c(control="#69b3a2", exclosure="grey"))+
  geom_vline(xintercept=0, lty=2)+ggtitle("PB (h=12)")+facet_wrap(~plot)

pb_h12

ggarrange(pb_h1,pb_h6,pb_h12, common.legend = T)

#alternative: density plots####

pb_h1=ggdensity(pbevals1, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-100, 40)

pb_h6=ggdensity(pbevals6, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), rug=T, add="mean", xlab=F,size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-100, 40)

pb_h12=ggdensity(pbevals12, x="score_diff", color="plot",fill="plot", size=1,
                 palette=c("#69b3a2", "grey"), rug=T, add="mean", xlab="RMSE difference (matched-mismatched)",
                 main="h=12")+geom_vline(xintercept=0)+xlim(-100, 40)

pbh=ggarrange(pb_h1,pb_h6,pb_h12, common.legend = T, nrow=3)
annotate_figure(pbh, top = text_grob("C. baileyi", face = "bold", size = 14))

###RMSE~moon#####

#control
colors <- c("same" = "darkblue", "switched" = "darkred")

#h=1#
pbevals1_control=pbevals1%>%filter(plot=="control")
pbevals6_control=pbevals6%>%filter(plot=="control")
pbevals12_control=pbevals12%>%filter(plot=="control")

pbevals1_exclosure=pbevals1%>%filter(plot=="exclosure")
pbevals6_exclosure=pbevals6%>%filter(plot=="exclosure")
pbevals12_exclosure=pbevals12%>%filter(plot=="exclosure")

pbevals_moon1_cont=ggplot()+
  geom_line(data=pbevals1_control, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=pbevals1_control, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PB control (h=1)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

pbevals_moon1_cont

pbevals_moon6_cont=ggplot()+
  geom_line(data=pbevals6_control, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=pbevals6_control, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PB control (h=6)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

pbevals_moon6_cont

pbevals_moon12_cont=ggplot()+
  geom_line(data=pbevals12_control, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=pbevals12_control, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PB control (h=12)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

pbevals_moon12_cont

ggarrange(pbevals_moon1_cont, pbevals_moon6_cont, pbevals_moon12_cont, common.legend = T, nrow=1, ncol=3)

#exclosure
pbevals_moon1_excl=ggplot()+
  geom_line(data=pbevals1_exclosure, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=pbevals1_exclosure, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PB exclosure (h=1)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

pbevals_moon1_excl

pbevals_moon6_excl=ggplot()+
  geom_line(data=pbevals6_exclosure, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=pbevals6_exclosure, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PB exclosure (h=6)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

pbevals_moon6_excl

pbevals_moon12_excl=ggplot()+
  geom_line(data=pbevals12_exclosure, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=pbevals12_exclosure, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PB exclosure (h=12)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

pbevals_moon12_excl

ggarrange(pbevals_moon1_excl, pbevals_moon6_excl, pbevals_moon12_excl, common.legend = T, nrow=1, ncol=3)