#Fig. 1: time series####

#DM data

dmcont_dat=rodent_data%>%
  mutate(abundance=rowSums(.[6:8]))%>%
  select(newmoonnumber,treatment, abundance)%>%
  replace_na(list(treatment='control'))%>%
  filter(treatment=="control")

dmexcl_dat=rodent_data%>%
  mutate(abundance=rowSums(.[6:8]))%>%
  select(newmoonnumber,treatment, abundance)%>%
  replace_na(list(treatment='exclosure'))%>%
  filter(treatment=="exclosure")

dmcont_dat1=dmcont_dat%>%filter(!newmoonnumber<278, !newmoonnumber>526)%>%mutate(species="Dipodomys spp.")
dmexcl_dat1=dmexcl_dat%>%filter(!newmoonnumber<278, !newmoonnumber>526)%>%mutate(species="Dipodomys spp.")

dmcont_dat1$abundance=round_na.interp(dmcont_dat1$abundance)
dmexcl_dat1$abundance=round_na.interp(dmexcl_dat1$abundance)

dmdat=rbind(dmcont_dat1, dmexcl_dat1)

#PB data
pbcont_dat=rodent_data%>%
  select(newmoonnumber,treatment, PB)%>%rename("abundance"="PB")%>%
  replace_na(list(treatment='control'))%>%
  filter(treatment=="control")

pbexcl_dat=rodent_data%>%
  select(newmoonnumber,treatment, PB)%>%rename("abundance"="PB")%>%
  replace_na(list(treatment='exclosure'))%>%
  filter(treatment=="exclosure")

pbcont_dat1=pbcont_dat%>%mutate(species="C. baileyi")%>%
  filter(!newmoonnumber<278, !newmoonnumber>526)

pbexcl_dat1=pbexcl_dat%>%mutate(species="C. baileyi")%>%
  filter(!newmoonnumber<278, !newmoonnumber>526)

pbcont_dat1$abundance=round_na.interp(pbcont_dat1$abundance)
pbexcl_dat1$abundance=round_na.interp(pbexcl_dat1$abundance)

pbdat=rbind(pbcont_dat1, pbexcl_dat1)

#PP data
ppcont_dat=rodent_data%>%
  select(newmoonnumber,treatment, PP)%>%rename("abundance"="PP")%>%
  replace_na(list(treatment='control'))%>%
  filter(treatment=="control")

ppexcl_dat=rodent_data%>%
  select(newmoonnumber,treatment, PP)%>%rename("abundance"="PP")%>%
  replace_na(list(treatment='exclosure'))%>%
  filter(treatment=="exclosure")

ppcont_dat1=ppcont_dat%>%mutate(species="C. penicillatus")%>%
  filter(!newmoonnumber<278, !newmoonnumber>526)

ppexcl_dat1=ppexcl_dat%>%mutate(species="C. penicillatus")%>%
  filter(!newmoonnumber<278, !newmoonnumber>526)

ppcont_dat1$abundance=round_na.interp(ppcont_dat1$abundance)
ppexcl_dat1$abundance=round_na.interp(ppexcl_dat1$abundance)

ppdat=rbind(ppcont_dat1, ppexcl_dat1)

spdat=rbind(dmdat, pbdat, ppdat)

spdat1=spdat%>%mutate(trt=case_when(treatment=="control"~"control",
                                    treatment=="exclosure"~"removal"))


#plot
pp11=ggplot(spdat1)+ geom_line(mapping=aes(x=newmoonnumber, y=abundance, col=species))+
  theme(legend.text = element_text(face="italic"),  
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  facet_wrap(~trt, nrow=2)+
  xlab("")+
  geom_bracket(
    xmin = 280 , xmax = 396, y.position = 60,
    label = paste("italic(C.baileyi)"), type="expression"
  )+
  geom_bracket(
    xmin = 411 , xmax = 526, y.position = 60,
    label = paste("italic(C.penicillatus)"), type="expression"
  )

#PP results####

##Fig. 2: coefficient overlap####
int_ov=ggdensity(coef_df_PP, x="intercept", color="treatment",fill="treatment", 
                 palette=c("#69b3a2", "grey"), rug=T, ylab=F, add="median",xlab=F, size=1
)+geom_vline(xintercept=0)+
  scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  
  theme(aspect.ratio=1,axis.text.x = element_text(angle=90), axis.text.y = element_blank())

ar1_ov=ggdensity(coef_df_PP, x="beta1", color="treatment",fill="treatment", 
                 palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+ 
  theme(aspect.ratio=1,axis.text.x = element_text(angle=90), axis.text.y = element_blank())

ar12_ov=ggdensity(coef_df_PP, x="beta12", color="treatment",fill="treatment",
                  palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  
  theme(aspect.ratio=1,axis.text.x = element_text(angle=90), axis.text.y = element_blank())

temp_ov=ggdensity(coef_df_PP, x="temp", color="treatment",fill="treatment", 
                  palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(aspect.ratio=1,axis.text.x = element_text(angle=90), axis.text.y = element_blank())

wprec_ov=ggdensity(coef_df_PP, x="warm_precip", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  scale_y_continuous(n.breaks=3)+
  theme(aspect.ratio=1,axis.text.x = element_text(angle=90), axis.text.y = element_blank())

cprec_ov=ggdensity(coef_df_PP, x="cool_precip", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  
  theme(aspect.ratio=1,axis.text.x = element_text(angle=90), axis.text.y = element_blank())

coef_ov=ggarrange(int_ov, ar1_ov, ar12_ov, temp_ov, wprec_ov, cprec_ov, common.legend = T, nrow=1, legend="bottom")+
  theme(plot.margin = margin(0,0,0,0, "cm"))+annotate("text", label="C. penicillatus", fontface="italic", size=14)

coef_ov1=coef_ov%>%gridExtra::grid.arrange(get_legend(coef_ov), heights = unit(c(70, 5), "mm"))

annotate_figure(coef_ov1, left = text_grob("C. penicillatus", face = "italic", size = 16, rot=90))

##Fig. 3: forecasts####

preds_plot1=ggplot()+
  geom_line(data=ppcont_covs_dat, aes(y=abundance, x=newmoonnumber), size=0.75)+
  geom_line(data=PPpreds_cont_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4)+
  geom_line(data=PPpreds_cont_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(colour="configuration")+
  ggtitle("control")+xlab("sampling number (new moon number)")+xlab("")+
  scale_colour_manual(values=c(same="darkblue", switched="darkred"), labels=c("non-transferred", "transferred"))

preds_plot1
annotate_figure(preds_plot1,left = text_grob("C. penicillatus", face = "italic", size = 16, rot=90))

preds_plot2=ggplot()+geom_line(data=ppexcl_covs_dat, aes(y=abundance, x=newmoonnumber), size=0.75)+
  geom_line(data=PPpreds_excl_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4)+
  geom_line(data=PPpreds_excl_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4)+
  theme_classic()+labs(colour="configuration")+xlab("sampling number (new moon number)")+
  ggtitle("removal")+
  scale_colour_manual(values=c(same="darkblue", switched="darkred"), labels=c("non-transferred", "transferred"))

preds_plot2
annotate_figure(preds_plot2,left = text_grob("C. penicillatus", face = "italic", size = 16, rot=90))

pp_predsplot=ggarrange(preds_plot1, preds_plot2, common.legend = T, ncol=1, nrow=2, legend="bottom")
annotate_figure(pp_predsplot,left = text_grob("C. penicillatus", face = "italic", size = 16, rot=90))

##Fig. 4: RMSE####
pp_h1=ggdensity(ppevals1, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), ylab=F,rug=T, add="median",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-25, 25)+
  annotate("text", x=-15, y=0.40, label="non-transferred better fit")+
  annotate("text", x=13, y=0.40, label="transferred better fit")+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(axis.text.y = element_blank())

pp_h6=ggdensity(ppevals6, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), rug=T,ylab=F, add="median", xlab=F,size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-25, 25)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(axis.text.y = element_blank())

pp_h12=ggdensity(ppevals12, x="score_diff", color="plot",fill="plot", size=1,
                 palette=c("#69b3a2", "grey"),ylab=F, rug=T, add="median", xlab="RMSE difference (non-transferred - transferred)",
                 main="h=12")+geom_vline(xintercept=0)+xlim(-25, 25)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(axis.text.y = element_blank())

pph=ggarrange(pp_h1,pp_h6,pp_h12, common.legend = T, nrow=3, legend="bottom")
annotate_figure(pph, top = text_grob("C. penicillatus", face = "italic", size = 16), 
                left= text_grob("RMSE", size=16, rot=90))

#Brier scores###

ppbr1=ggdensity(pp_briers1, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"),ylab=F, rug=T, add="median",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-0.20, .20)+
  annotate("text", x=-.12, y=70, label="non-transferred better fit")+
  annotate("text", x=0.10, y=70, label="transferred better fit")+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(axis.text.y = element_blank())

ppbr2=ggdensity(pp_briers6, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), ylab=F, rug=T, add="median",xlab=F, size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-0.20, .20)+
  theme(axis.text.y = element_blank())+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

ppbr3=ggdensity(pp_briers12, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"),ylab=F, rug=T, add="median", size=1,
                main="h=12", xlab="Brier score difference (non-transferred - transferred)")+
  geom_vline(xintercept=0)+xlim(-0.2, 0.2)+
  theme(axis.text.y = element_blank())+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

ppbh=ggarrange(ppbr1,ppbr2,ppbr3, common.legend = T, nrow=3, legend="bottom")
annotate_figure(ppbh, top = text_grob("C. penicillatus", face = "italic", size = 16), 
                left= text_grob("Brier score", size=16, rot=90))

#PB results####

##Fig. 2: coefficient overlap####
pbint_ov=ggdensity(coef_df_PB, x="intercept", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="median",xlab="intercept", size=1
)+geom_vline(xintercept=0)+ylab("")+
  scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(aspect.ratio = 1,axis.text.x = element_text(angle=90), axis.text.y = element_blank())

pbar1_ov=ggdensity(coef_df_PB, x="beta1", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="median",xlab="AR (1)", size=1, ylab=F
)+geom_vline(xintercept=0)+
  scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(aspect.ratio = 1, axis.text.x = element_text(angle=90), axis.text.y = element_blank())

pbar12_ov=ggdensity(coef_df_PB, x="beta12", color="treatment",fill="treatment",
                    palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=" AR (12)", size=1, ylab=F
)+geom_vline(xintercept=0)+
  scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(aspect.ratio = 1, axis.text.x = element_text(angle=90), axis.text.y = element_blank())

pbtemp_ov=ggdensity(coef_df_PB, x="temp", color="treatment",fill="treatment", 
                    palette=c("#69b3a2", "grey"), rug=T, add="median",xlab="mean temperature (lag=1)", size=1, ylab=F
)+geom_vline(xintercept=0)+
  scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(aspect.ratio = 1,axis.text.x = element_text(angle=90), axis.text.y = element_blank())

pbwprec_ov=ggdensity(coef_df_PB, x="warm_precip", color="treatment",fill="treatment", 
                     palette=c("#69b3a2", "grey"), rug=T, add="median",xlab="warm precipitation", size=1, ylab=F
)+geom_vline(xintercept=0)+
  scale_y_continuous(n.breaks=3)+
  scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(aspect.ratio = 1,axis.text.x = element_text(angle=90), axis.text.y = element_blank())

pbcprec_ov=ggdensity(coef_df_PB, x="cool_precip", color="treatment",fill="treatment", 
                     palette=c("#69b3a2", "grey"), rug=T, add="median",xlab="cool precipitation", size=1, ylab=F
)+geom_vline(xintercept=0)+
  scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(aspect.ratio = 1, axis.text.x = element_text(angle=90), axis.text.y = element_blank())

pbcoef_ov=ggarrange(pbint_ov, pbar1_ov, pbar12_ov, pbtemp_ov, pbwprec_ov, pbcprec_ov, common.legend = T, nrow=1, legend="bottom")+
  theme(plot.margin = margin(0,0,0,0, "cm"))+annotate("text", label="C. baileyi", fontface="italic", size=14)

pbcoef_ov1=pbcoef_ov%>%gridExtra::grid.arrange(get_legend(pbcoef_ov), heights = unit(c(70, 5), "mm"))

annotate_figure(pbcoef_ov1, left = text_grob("C. baileyi", face = "italic", size = 16, rot=90))

##Fig. 3: forecasts####

pbpreds_plot1=ggplot()+
  geom_line(data=pbcont_covs_dat, aes(y=abundance, x=newmoonnumber), size=0.75)+
  geom_line(data=PBpreds_cont_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4)+
  geom_line(data=PBpreds_cont_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  labs(colour="configuration")+
  ggtitle("control")+xlab("sampling number (new moon number)")+xlab("")+
  scale_colour_manual(values=c(same="darkblue", switched="darkred"), labels=c("non-transferred", "transferred"))

pbpreds_plot1
annotate_figure(pbpreds_plot1,left = text_grob("C. baileyi", face = "italic", size = 16, rot=90))

pbpreds_plot2=ggplot()+geom_line(data=pbexcl_covs_dat, aes(y=abundance, x=newmoonnumber), size=0.75)+
  geom_line(data=PBpreds_excl_same, aes(y=preds_same, x=moon, group=m,color="same"), alpha=0.4)+
  geom_line(data=PBpreds_excl_switch, aes(y=preds_switch, x=moon, group=m, color="switched"), alpha=0.4)+
  theme_classic()+labs(colour="configuration")+xlab("sampling number (new moon number)")+
  ggtitle("removal")+
  scale_colour_manual(values=c(same="darkblue", switched="darkred"), labels=c("non-transferred", "transferred"))

pbpreds_plot2
annotate_figure(pbpreds_plot2,left = text_grob("C. baileyi", face = "italic", size = 16, rot=90))

pb_predsplot=ggarrange(pbpreds_plot1, pbpreds_plot2, common.legend = T, ncol=1, nrow=2, legend="bottom")
annotate_figure(pb_predsplot,left = text_grob("C. baileyi", face = "italic", size = 16, rot=90))

##Fig. 4: RMSE####
pb_h1=ggdensity(pbevals1, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), rug=T, ylab=F,add="median",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-45, 45)+
  annotate("text", x=-28, y=0.12, label="non-transferred better fit")+
  annotate("text", x=26, y=0.12, label="transferred better fit")+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(axis.text.y = element_blank())

pb_h6=ggdensity(pbevals6, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), rug=T, ylab=F, add="median", xlab=F,size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-45, 45)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  theme(axis.text.y = element_blank())

pb_h12=ggdensity(pbevals12, x="score_diff", color="plot",fill="plot", size=1,ylab=F,
                 palette=c("#69b3a2", "grey"), rug=T, add="median", xlab="RMSE difference (non-transferred - transferred)",
                 main="h=12")+geom_vline(xintercept=0)+xlim(-45, 45)+
  theme(axis.text.y = element_blank())+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

pbhh=ggarrange(pb_h1,pb_h6,pb_h12, common.legend = T, nrow=3, legend="bottom")
annotate_figure(pbhh, top = text_grob("C. baileyi", face = "italic", size = 16), 
                left= text_grob("RMSE", size=16, rot=90))

#Brier scores###

pbbr1=ggdensity(pb_briers1, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T,ylab=F, add="median",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-0.4, 0.4)+
  annotate("text", x=-0.25, y=13, label="non-transferred better fit")+
  annotate("text", x=.25, y=13, label="transferred better fit")+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+theme(axis.text.y = element_blank())

pbbr2=ggdensity(pb_briers6, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F,ylab=F, size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-0.4, 0.4)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+theme(axis.text.y = element_blank())

pbbr3=ggdensity(pb_briers12, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T,ylab=F, add="median", size=1,
                main="h=12", xlab="Brier score difference (non-transferred - transferred)")+
  geom_vline(xintercept=0)+xlim(-0.4, 0.4)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+theme(axis.text.y = element_blank())

pbbh=ggarrange(pbbr1,pbbr2,pbbr3, common.legend = T, nrow=3, legend="bottom")
annotate_figure(pbbh, top = text_grob("C. baileyi", face = "italic", size = 16), 
                left= text_grob("Brier score", size=16, rot=90))
