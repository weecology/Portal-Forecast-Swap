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

dmcont_covs_dat=dmcont_dat%>%filter(!newmoonnumber<278, !newmoonnumber>526)%>%mutate(species="Dipodomys spp.")
dmexcl_covs_dat=dmexcl_dat%>%filter(!newmoonnumber<278, !newmoonnumber>526)%>%mutate(species="Dipodomys spp.")

dmcont_covs_dat$abundance=round_na.interp(dmcont_covs_dat$abundance)
dmexcl_covs_dat$abundance=round_na.interp(dmexcl_covs_dat$abundance)

dmdat=rbind(dmcont_covs_dat, dmexcl_covs_dat)

#PB data
pbcont_covs_dat=pbcont_covs_dat%>%mutate(species="C. baileyi", treatment="control")%>%select(newmoonnumber, treatment, abundance, species)
pbexcl_covs_dat=pbexcl_covs_dat%>%mutate(species="C. baileyi", treatment="exclosure")%>%select(newmoonnumber, treatment, abundance, species)

pbdat=rbind(pbcont_covs_dat, pbexcl_covs_dat)

#PP data
ppcont_covs_dat=ppcont_covs_dat%>%mutate(species="C. penicillatus", treatment="control")%>%select(newmoonnumber, treatment, abundance, species)
ppexcl_covs_dat=ppexcl_covs_dat%>%mutate(species="C. penicillatus", treatment="exclosure")%>%select(newmoonnumber, treatment, abundance, species)

ppdat=rbind(ppcont_covs_dat, ppexcl_covs_dat)

spdat=rbind(dmdat, pbdat, ppdat)

spdat_cont=spdat%>%filter(treatment=="control")
spdat_excl=spdat%>%filter(treatment=="exclosure")

#plot
pp11=ggplot()+theme(legend.text=element_text(face="italic"),
                     panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line(data=spdat_cont, aes(x=newmoonnumber, y=abundance, col=species))+
  ggtitle("control")+
  xlab("")

pp12=ggplot()+theme(legend.text=element_text(face="italic"),
                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                    panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  geom_line(data=spdat_excl, aes(x=newmoonnumber, y=abundance, col=species))+
  xlab("sampling number (new moon number)")+ggtitle("removal")

ggarrange(pp11, pp12, common.legend = T, legend = "bottom", nrow=2)

#PP results####

##Fig. 2: coefficient overlap####
int_ov=ggdensity(coef_df_PP, x="intercept", color="treatment",fill="treatment", 
                 palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

ar1_ov=ggdensity(coef_df_PP, x="beta1", color="treatment",fill="treatment", 
                 palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

ar12_ov=ggdensity(coef_df_PP, x="beta12", color="treatment",fill="treatment",
                  palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

temp_ov=ggdensity(coef_df_PP, x="temp", color="treatment",fill="treatment", 
                  palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

wprec_ov=ggdensity(coef_df_PP, x="warm_precip", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  scale_y_continuous(n.breaks=3)+
  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

cprec_ov=ggdensity(coef_df_PP, x="cool_precip", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

coef_ov=ggarrange(int_ov, ar1_ov, ar12_ov, temp_ov, wprec_ov, cprec_ov, common.legend = T, nrow=1, legend="bottom")+
  theme(plot.margin = margin(0,0,0,0, "cm"))+annotate("text", label="C. penicillatus", fontface="italic", size=14)

annotate_figure(coef_ov, left = text_grob("C. penicillatus", face = "italic", size = 16, rot=90))

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
                palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-25, 25)+
  annotate("segment", x=-5, y= 0.45, xend=-23, yend=0.45, col="black", size=1, arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=5, y= 0.45, xend=23, yend=0.45, col="black", size=1, arrow=arrow(length=unit(0.3, "cm")))+
  annotate("text", x=-15, y=0.40, label="non-transferred better fit")+
  annotate("text", x=13, y=0.40, label="transferred better fit")+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

pp_h6=ggdensity(ppevals6, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), rug=T, add="median", xlab=F,size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-25, 25)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

pp_h12=ggdensity(ppevals12, x="score_diff", color="plot",fill="plot", size=1,
                 palette=c("#69b3a2", "grey"), rug=T, add="median", xlab="RMSE difference (non-transferred - transferred)",
                 main="h=12")+geom_vline(xintercept=0)+xlim(-25, 25)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

pph=ggarrange(pp_h1,pp_h6,pp_h12, common.legend = T, nrow=3, legend="bottom")
annotate_figure(pph,left = text_grob("C. penicillatus", face = "italic", size = 16, rot=90))

##Fig. 5: Brier scores####

ppbr1=ggdensity(pp_briers1, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-0.20, .20)+
  annotate("segment", x=-0.05, y= 98, xend=-0.18, yend=98, col="black", size=1, arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", xend=0.18, y= 98, x=.05, yend=98, col="black", size=1, arrow=arrow(length=unit(0.3, "cm")))+
  annotate("text", x=-.12, y=87, label="non-transferred better fit")+
  annotate("text", x=0.10, y=87, label="transferred better fit")+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

ppbr2=ggdensity(pp_briers6, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-0.20, .20)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

ppbr3=ggdensity(pp_briers12, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T, add="median", size=1,
                main="h=12", xlab="Brier score difference (non-transferred - transferred)")+
  geom_vline(xintercept=0)+xlim(-0.2, 0.2)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

ppbh=ggarrange(ppbr1,ppbr2,ppbr3, common.legend = T, nrow=3, legend="bottom")
annotate_figure(ppbh,left = text_grob("C. penicillatus", face = "italic", size = 16, rot=90))

#PB results####

##Fig. 2: coefficient overlap####
pbint_ov=ggdensity(coef_df_PB, x="intercept", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="median",xlab="intercept", size=1
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

pbar1_ov=ggdensity(coef_df_PB, x="beta1", color="treatment",fill="treatment", 
                   palette=c("#69b3a2", "grey"), rug=T, add="median",xlab="AR (1)", size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

pbar12_ov=ggdensity(coef_df_PB, x="beta12", color="treatment",fill="treatment",
                    palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=" AR (12)", size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

pbtemp_ov=ggdensity(coef_df_PB, x="temp", color="treatment",fill="treatment", 
                    palette=c("#69b3a2", "grey"), rug=T, add="median",xlab="mean temperature (lag=1)", size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

pbwprec_ov=ggdensity(coef_df_PB, x="warm_precip", color="treatment",fill="treatment", 
                     palette=c("#69b3a2", "grey"), rug=T, add="median",xlab="warm precipitation", size=1, ylab=F
)+geom_vline(xintercept=0)+scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+
  scale_y_continuous(n.breaks=3)+
  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

pbcprec_ov=ggdensity(coef_df_PB, x="cool_precip", color="treatment",fill="treatment", 
                     palette=c("#69b3a2", "grey"), rug=T, add="median",xlab="cool precipitation", size=1, ylab=F
)+geom_vline(xintercept=0)+
  scale_x_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))+  theme(axis.text.x = element_text(angle=90), axis.text.y = element_text(angle=90))

pbcoef_ov=ggarrange(pbint_ov, pbar1_ov, pbar12_ov, pbtemp_ov, pbwprec_ov, pbcprec_ov, common.legend = T, nrow=1, legend="bottom")+
  theme(plot.margin = margin(0,0,0,0, "cm"))+annotate("text", label="C. baileyi", fontface="italic", size=14)

annotate_figure(pbcoef_ov, left = text_grob("C. baileyi", face = "italic", size = 16, rot=90))

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
                palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-50, 50)+
  annotate("segment", x=-15, y= 0.15, xend=-45, yend=0.15, col="black", size=1, arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=15, y= 0.15, xend=45, yend=0.15, col="black", size=1, arrow=arrow(length=unit(0.3, "cm")))+
  annotate("text", x=-28, y=0.14, label="non-transferred better fit")+
  annotate("text", x=26, y=0.14, label="transferred better fit")+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

pb_h6=ggdensity(pbevals6, x="score_diff", color="plot",fill="plot", 
                palette=c("#69b3a2", "grey"), rug=T, add="median", xlab=F,size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-50, 50)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

pb_h12=ggdensity(pbevals12, x="score_diff", color="plot",fill="plot", size=1,
                 palette=c("#69b3a2", "grey"), rug=T, add="median", xlab="RMSE difference (non-transferred - transferred)",
                 main="h=12")+geom_vline(xintercept=0)+xlim(-50, 50)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

pbhh=ggarrange(pb_h1,pb_h6,pb_h12, common.legend = T, nrow=3, legend="bottom")
annotate_figure(pbhh,left = text_grob("C. baileyi", face = "italic", size = 16, rot=90))

##Fig. 5: Brier scores####

pbbr1=ggdensity(pb_briers1, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1,
                main="h=1")+geom_vline(xintercept=0)+xlim(-0.4, 0.4)+
  annotate("segment", x=-.15, y= 14, xend=-.35, yend=14, col="black", size=1, arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=.15, y= 14, xend=.35, yend=14, col="black", size=1, arrow=arrow(length=unit(0.3, "cm")))+
  annotate("text", x=-0.25, y=13, label="matched better fit")+
  annotate("text", x=.25, y=13, label="mismatched better fit")+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

pbbr2=ggdensity(pb_briers6, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T, add="median",xlab=F, size=1,
                main="h=6")+geom_vline(xintercept=0)+xlim(-0.4, 0.4)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

pbbr3=ggdensity(pb_briers12, x="brier_diff", color="treatment",fill="treatment", 
                palette=c("#69b3a2", "grey"), rug=T, add="median", size=1,
                main="h=12", xlab="Brier score difference (non-transferred - transferred)")+
  geom_vline(xintercept=0)+xlim(-0.4, 0.4)+
  scale_y_continuous(n.breaks=3, labels = scales::label_number(accuracy = 0.01))

pbbh=ggarrange(pbbr1,pbbr2,pbbr3, common.legend = T, nrow=3, legend="bottom")
annotate_figure(pbbh,left = text_grob("C. baileyi", face = "italic", size = 16, rot=90))
