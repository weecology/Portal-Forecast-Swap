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

#coefficient overlap####
int_ov=ggdensity(coef_df_PB, x="intercept", color="treatment",fill="treatment", 
                 palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                 main="intercept")+geom_vline(xintercept=0)

ar1_ov=ggdensity(coef_df_PB, x="beta1", color="treatment",fill="treatment", 
                 palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                 main="AR (1)")+geom_vline(xintercept=0)

ar13_ov=ggdensity(coef_df_PB, x="beta13", color="treatment",fill="treatment", 
                  palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                  main="AR(13)")+geom_vline(xintercept=0)

temp_ov=ggdensity(coef_df_PB, x="temp", color="treatment",fill="treatment", 
                  palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                  main="mean temperature (lag=1)")+geom_vline(xintercept=0)
wprec_ov=ggdensity(coef_df_PB, x="warm_precip", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                   main="warm precipitation")+geom_vline(xintercept=0)
cprec_ov=ggdensity(coef_df_PB, x="cool_precip", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                   main="cool precipitation")+geom_vline(xintercept=0)

coef_ov=ggarrange(int_ov, ar1_ov, ar13_ov, temp_ov, wprec_ov, cprec_ov, common.legend = T)
annotate_figure(coef_ov, top=text_grob("C. baileyi", face="bold", size=14))

##forecasts####

#h=12
#control-control
pbcont_preds_same=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$preds_same), get_dat_same)

#control dat-exclosure mod
pbcont_preds_switch=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$preds_switch), get_dat_switch)

m=rep(seq(1:59), each=13)
code1="same"
code2="switched"

PBpreds_cont_same=do.call(rbind.data.frame, pbcont_preds_same)
PBpreds_cont_same=cbind(PBpreds_cont_same, m, code1)
PBpreds_cont_switch=do.call(rbind.data.frame, pbcont_preds_switch)
PBpreds_cont_switch=cbind(PBpreds_cont_switch, m, code2)

pbpreds_plot1=ggplot()+
  geom_line(data=pbcont_covs_ro, aes(y=abundance, x=newmoonnumber), size=0.75)+
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

pbpreds_plot2=ggplot()+geom_line(data=pbexcl_covs_ro, aes(y=abundance, x=newmoonnumber), size=0.75)+
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

#h=1
pbcont_evals1_diff=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$evals_same, PBcontrol_dat$evals_switch ,PBcontrol_dat$id), get_evals1_diff)
pbevals1_cont_diff=do.call(rbind.data.frame, pbcont_evals1_diff)%>%mutate(plot="control")

pbexcl_evals1_diff=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$evals_same, PBexclosure_dat$evals_switch ,PBexclosure_dat$id), get_evals1_diff)
pbevals1_excl_diff=do.call(rbind.data.frame, pbexcl_evals1_diff)%>%mutate(plot="exclosure")

pbevals1=rbind(pbevals1_cont_diff, pbevals1_excl_diff)
pbevals1$newmoon=as.integer(pbevals1$newmoon)
pbevals1$h=as.integer(pbevals1$h)
pbevals1$score_same=as.numeric(pbevals1$score_same)
pbevals1$score_switch=as.numeric(pbevals1$score_switch)
pbevals1$score_diff=as.numeric(pbevals1$score_diff)

str(pbevals1)

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

#h=13
pbcont_evals13_diff=pmap(list(PBcontrol_dat$splits,PBcontrol_dat$evals_same, PBcontrol_dat$evals_switch ,PBcontrol_dat$id), get_evals13_diff)
pbevals13_cont_diff=do.call(rbind.data.frame, pbcont_evals13_diff)%>%mutate(plot="control")

pbexcl_evals13_diff=pmap(list(PBexclosure_dat$splits,PBexclosure_dat$evals_same, PBexclosure_dat$evals_switch ,PBexclosure_dat$id), get_evals13_diff)
pbevals13_excl_diff=do.call(rbind.data.frame, pbexcl_evals13_diff)%>%mutate(plot="exclosure")

pbevals13=rbind(pbevals13_cont_diff, pbevals13_excl_diff)
pbevals13$newmoon=as.integer(pbevals13$newmoon)
pbevals13$h=as.integer(pbevals13$h)
pbevals13$score_same=as.numeric(pbevals13$score_same)
pbevals13$score_switch=as.numeric(pbevals13$score_switch)
pbevals13$score_diff=as.numeric(pbevals13$score_diff)

str(pbevals13)

#plot
pb_h1=ggdensity(pbevals1, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-250, 250)

pb_h6=ggdensity(pbevals6, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), rug=T, add="mean", xlab=F,size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-250, 250)

pb_h13=ggdensity(pbevals13, x="score_diff", color="plot",fill="plot", size=1,
                 palette=c("#69b3a2", "grey"), rug=T, add="mean", xlab="RMSE difference (matched-mismatched)",
                 main="h=13")+geom_vline(xintercept=0)+xlim(-250, 250)

pbh=ggarrange(pb_h1,pb_h6,pb_h13, common.legend = T, nrow=3)
annotate_figure(pbh, top = text_grob("C. baileyi", face = "bold", size = 14))

###Brier score####

#control###

#combine predictions on same and switched models for control data
pb_preds_control=left_join(PBpreds_cont_same, PBpreds_cont_switch, by=c("moon", "holdout", "m"))

#calculate brier score for same models
pb_brier_cont1=scoring(pred=pb_preds_control$preds_same, response=pb_preds_control$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine same model predictions and brier score dataframes
pb_brier_cont_same=cbind(pb_preds_control, pb_brier_cont1)%>%rename(quadratic_same=quadratic)

#calculate brier score for switched models
pb_brier_cont2=scoring(pred=pb_preds_control$preds_switch, response=pb_preds_control$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine switched model predictions and brier score dataframes
pb_brier_cont_switch=cbind(pb_preds_control, pb_brier_cont2)%>%rename(quadratic_switch=quadratic)

#calculate brier score differences
pb_brier_control=left_join(pb_brier_cont_same, pb_brier_cont_switch)%>%mutate(treatment="control",brier_diff=quadratic_same-quadratic_switch)

#select relevant rows with horizons 1,6 13
pb_brier_control1=pb_brier_control%>%group_by(m)%>%slice(which(row_number()%in%c(1,6,13)))

#subdivide into each horizon for easier plotting
pbb1=pb_brier_control1%>%group_by(m)%>%slice_head(n=1)%>%mutate(horizon="1")
pbb6=pb_brier_control1%>%group_by(m)%>%slice(which(row_number()==2))%>%mutate(horizon="6")
pbb13=pb_brier_control1%>%group_by(m)%>%slice(which(row_number()==3))%>%mutate(horizon="13")

#combine all 3 horizons
pbb=rbind(pbb1,pbb6,pbb13)

#exclosure###

#combine predictions on same and switched models for exclosure data
pb_preds_exclosure=left_join(PBpreds_excl_same, PBpreds_excl_switch, by=c("moon", "holdout", "m"))

#calculate brier score for same models
pb_brier_excl1=scoring(pred=pb_preds_exclosure$preds_same, response=pb_preds_exclosure$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine same model predictions and brier score dataframes
pb_brier_excl_same=cbind(pb_preds_exclosure, pb_brier_excl1)%>%rename(quadratic_same=quadratic)

#calculate brier score for switched models
pb_brier_excl2=scoring(pred=pb_preds_exclosure$preds_switch, response=pb_preds_exclosure$holdout,distr="nbinom", distrcoefs=2, individual=TRUE,
                       cutoff=1000)%>%select(quadratic)

#combine switched model predictions and brier score dataframes
pb_brier_excl_switch=cbind(pb_preds_exclosure, pb_brier_excl2)%>%rename(quadratic_switch=quadratic)

#calculate brier score differences
pb_brier_exclosure=left_join(pb_brier_excl_same, pb_brier_excl_switch)%>%mutate(treatment="exclosure",brier_diff=quadratic_same-quadratic_switch)

#select relevant rows with horizons 1,6 13
pb_brier_exclosure1=pb_brier_exclosure%>%group_by(m)%>%slice(which(row_number()%in%c(1,6,13)))

#subdivide into each horizon for easier plotting
pbbx1=pb_brier_exclosure1%>%group_by(m)%>%slice_head(n=1)%>%mutate(horizon="1")
pbbx6=pb_brier_exclosure1%>%group_by(m)%>%slice(which(row_number()==2))%>%mutate(horizon="6")
pbbx13=pb_brier_exclosure1%>%group_by(m)%>%slice(which(row_number()==3))%>%mutate(horizon="13")

#combine all 3 horizons
pbbx=rbind(pbbx1,pbbx6, pbbx13)

#combine control and exclosure and filter out per horizon for plotting
pb_briers1=rbind(pbb,pbbx)%>%filter(horizon==1)
pb_briers6=rbind(pbb,pbbx)%>%filter(horizon==6)
pb_briers13=rbind(pbb,pbbx)%>%filter(horizon==13)

pbbp1=ggdensity(pb_briers1, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-0.50, .50)

pbbp2=ggdensity(pb_briers6, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T, add="mean",xlab=F, size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-0.50, .50)

pbbp3=ggdensity(pb_briers13, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T, add="mean", size=1,
                main="h=13", xlab="Brier score difference (matched-mismatched)")+geom_vline(xintercept=0)+xlim(-0.50, .50)

pbbh=ggarrange(pbbp1,pbbp2,pbbp3, common.legend = T, nrow=3)
annotate_figure(pbbh, top = text_grob("C. baileyi", face = "bold", size = 14))

#extra stuff;unused plot###
###RMSE~moon#####

#control
colors <- c("same" = "darkblue", "switched" = "darkred")

#h=1#
pbevals1_control=pbevals1%>%filter(plot=="control")
pbevals6_control=pbevals6%>%filter(plot=="control")
pbevals13_control=pbevals13%>%filter(plot=="control")

pbevals1_exclosure=pbevals1%>%filter(plot=="exclosure")
pbevals6_exclosure=pbevals6%>%filter(plot=="exclosure")
pbevals13_exclosure=pbevals13%>%filter(plot=="exclosure")

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

pbevals_moon13_cont=ggplot()+
  geom_line(data=pbevals13_control, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=pbevals13_control, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PB control (h=13)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

pbevals_moon13_cont

ggarrange(pbevals_moon1_cont, pbevals_moon6_cont, pbevals_moon13_cont, common.legend = T, nrow=1, ncol=3)

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

pbevals_moon13_excl=ggplot()+
  geom_line(data=pbevals13_exclosure, aes(x=newmoon, y=score_same, color="same"), alpha=0.4)+
  geom_line(data=pbevals13_exclosure, aes(x=newmoon, y=score_switch, color="switched"), alpha=0.4)+
  theme_classic()+ylab("RMSE")+ggtitle("PB exclosure (h=13)")+labs(colour="configuration")+
  scale_color_manual(values=c(same="darkblue", switched="darkred"), labels=c("same", "switched"))

pbevals_moon13_excl

ggarrange(pbevals_moon1_excl, pbevals_moon6_excl, pbevals_moon13_excl, common.legend = T, nrow=1, ncol=3)