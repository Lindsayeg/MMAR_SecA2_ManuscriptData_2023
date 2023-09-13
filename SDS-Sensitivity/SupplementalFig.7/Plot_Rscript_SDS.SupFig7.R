
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(ggplot2,tidyverse,plyr,dplyr,DescTools,multcomp,ggpubr,rstatix,ggprism,magrittr,patchwork,ggtext,tibble,ggforce,patchwork)

setwd("~/Desktop/")
data=data.frame(read.csv("2023.07_SDS-Sensitivity-Assays.csv",header=T))

#Make the strain levels a factor and order as desired
data$Strain <- factor(data$Strain, levels = c("WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT"))



########Bar plot with data points as dots#######

#Plot IFNB concentration by strain

L <- ggplot(data, aes(x=Strain, y=RelGrowth, fill=Strain)) + geom_point(size=0.1) + theme_classic() + theme(axis.text=element_text(size=10),
                                                                                   axis.title.y=element_text(size=10),
                                                                                   axis.title.x=element_blank(),
                                                                                   axis.text.x=element_text(angle=45, hjust = 1))
#Add data points as dots and bar graph for average values
#AveL <- L + stat_summary(geom = "bar", fun = "mean", fill=c(rep("Black",5),"#990066"), alpha=c(0.2)) + ylab("IFN-β Concentration (pg/ml)") 
P1 <- L + stat_summary(geom = "bar", fun = "mean", color=c("#990000", "#FFCC33", "#CC9900","#996600"), linewidth=0.4, fill=c("#990000", "#FFCC33", "#CC9900","#996600"), alpha=c(0.4)) +
  theme(axis.text.x=element_markdown(angle=45, hjust = 1), axis.title.x=element_blank(), axis.title.y=element_blank(), legend.title = element_blank(), legend.position="none") +
  scale_x_discrete(labels = c("Uninfected"="Uninfected", "WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*", "ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*")) +
  #scale_y_continuous(trans='log10') +
  coord_cartesian(ylim = c(0.0001, 0.035))
  #+facet_zoom(ylim = c(0, 0.1))

P1 <- P1 + geom_signif(annotation=c("***", "***", "***"),
                      y_position = c(0.032,0.032,0.032),
                      xmin=c(1,3,4), xmax=c(1,3,4), textsize = 3, 
                      tip_length = 0, vjust=0, size=0)

P1

P2 <- L + stat_summary(geom = "bar", fun = "mean", color=c("#990000", "#FFCC33", "#CC9900","#996600"), linewidth=0.4, fill=c("#990000", "#FFCC33", "#CC9900","#996600"), alpha=c(0.4)) + ylab("Bacterial growth in SDS \n relative to no SDS") +
  theme(axis.text.x=element_markdown(angle=45, hjust = 1, size=12), axis.title.y=element_text(size=12), axis.title.x=element_blank(), legend.title = element_blank(), legend.position="none") +
  scale_x_discrete(labels = c("Uninfected"="Uninfected", "WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*", "ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*")) +
  #scale_y_continuous(trans='log10') +
  coord_cartesian(ylim = c(0.0001, 10))
#+facet_zoom(ylim = c(0, 0.1))

P2 + patchwork::inset_element(P1, 0.3, 0.3, 1, 1)



##CODE to Save below as a part of Statistics






#
##
###
####
#####
######
#######
########
#########
##########


####

###Use a Shapiro-Wilk normality test to determine if the data is normally distributed.
shapiro.test(data$RelGrowth)
shapiro.test(data$RelGrowth[data$Strain=="WT"])
shapiro.test(data$RelGrowth[data$Strain=="ΔsecA2"])
shapiro.test(data$RelGrowth[data$Strain=="ΔsecA2/secA2_MM"])
shapiro.test(data$RelGrowth[data$Strain=="ΔsecA2/secA2_MT"])



##### ALTERNATIVELY TRY A KRUSKAL-WALLIS test
kruskal.test(RelGrowth ~ Strain, data = data)

##As the p-value is less than 0.05, there is likely a significant interaction in the data set.
##Use a pair-wise comparison to determine which pair is significantly different.
pairwise.wilcox.test(data$RelGrowth, data$Strain, exact=F)


##Save plot
ggsave("./SDS-Sensitivity2.png", width = 5, height = 5)


























