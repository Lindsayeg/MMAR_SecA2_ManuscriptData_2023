
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(ggplot2,tidyverse,plyr,dplyr,DescTools,multcomp,ggpubr,rstatix,ggprism,magrittr,patchwork,ggtext,tibble)

setwd("~/Documents/ND_PhD/Writing Projects/Paper/MMAR_SecA2_ManuscriptData_2023/SDS-Sensitivity/")
data=data.frame(read.csv("SDS-Sensitivity_Counts.csv",header=T))

#Make the strain levels a factor and order as desired
data$Strain <- factor(data$Strain, levels = c("WT", "ΔsecA2", "ΔsecA2/secA2_MM"))



########Bar plot with data points as dots#######

#Plot IFNB concentration by strain

L <- ggplot(data, aes(x=Strain, y=RelGrowth, fill=Strain)) + geom_point(size=0.1) + theme_classic() + theme(axis.text=element_text(size=10),
                                                                                   axis.title.y=element_text(size=10),
                                                                                   axis.title.x=element_blank(),
                                                                                   axis.text.x=element_text(angle=45, hjust = 1))
#Add data points as dots and bar graph for average values
#AveL <- L + stat_summary(geom = "bar", fun = "mean", fill=c(rep("Black",5),"#990066"), alpha=c(0.2)) + ylab("IFN-β Concentration (pg/ml)") 
AveL <- L + stat_summary(geom = "bar", fun = "mean", color=c("#990000", "#FFCC33", "#CC9900"), linewidth=0.4, fill=c("#990000", "#FFCC33", "#CC9900"), alpha=c(0.4)) + ylab("Bacterial growth in SDS \n relative to no SDS") +
  theme(axis.text.x=element_markdown(angle=45, hjust = 1), axis.title.x=element_blank(), legend.title = element_blank(), legend.position="none") +
  scale_x_discrete(labels = c("Uninfected"="Uninfected", "WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*")) +
  #scale_y_continuous(trans='log10') +
  coord_cartesian(ylim = c(0.001, 6))


AveL

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
##Add One-Way ANOVA statistics
#Set Control - here I want everything compared to WT and not Uninfected so WT must be listed first
data$Strain <- factor(data$Strain, levels = c("WT","ΔsecA2", "ΔsecA2/secA2_MM"), ordered = TRUE)
#compute One-way ANVOA
res.aov <- aov(RelGrowth ~ Strain, data = data)

#Look at summary statistics
#Can see there is statistical significance but can't see between which variables
summary(res.aov)


#Run a Dunnett Test to determine which treatment are significant compared to the control (WT)
post_test <- glht(model = res.aov, linfct = mcp(Strain = "Dunnett"))

#Look at summary statistics
summary(post_test)


#Manually add ANOVA results (SEE BELOW)
p2 <- AveL + geom_signif(annotation=c("**", "ns"),
                         y_position = c(4,5),
                         xmin=c(1,1),
                         xmax=c(2,3), textsize = 4) 

p2



##Save plot
ggsave("./SDS-Sensitivity.png", width = 2.5, height = 2.5)




###Use a Shapiro-Wilk normality test to determine if the data is normally distributed.
shapiro.test(data$RelGrowth)
shapiro.test(data$RelGrowth[data$Strain=="WT"])
shapiro.test(data$RelGrowth[data$Strain=="ΔsecA2"])
shapiro.test(data$RelGrowth[data$Strain=="ΔsecA2/secA2_MM"])



##ALL Seem to be NORMAL (although the d_secA2 is very close!) so an ANVOA is okay?

##### ALTERNATIVELY TRY A KRUSKAL-WALLIS test
kruskal.test(RelGrowth ~ Strain, data = data)

##As the p-value is less than 0.05, there is likely a significant interaction in the data set.
##Use a pair-wise comparison to determine which pair is significantly different.
##If comparing to WT, this shows the d_secA2 strain is statistically significant from WT. The others are not, although d_esxBA is close.
pairwise.wilcox.test(data$RelGrowth, data$Strain, exact=F)









##Data and calculate SD, SE,and IC
DfDnew <- data %>% select(Strain, RelGrowth)
mysum <- DfDnew %>%
  group_by(Strain) %>% 
  dplyr::summarise(
    n=n(),
    mean= mean(RelGrowth), 
    sd = sd(RelGrowth)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))














