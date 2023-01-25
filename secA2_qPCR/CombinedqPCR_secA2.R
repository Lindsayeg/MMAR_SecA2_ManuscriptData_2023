## Install and load packages
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(ggplot2,tidyverse,plyr,dplyr,DescTools,multcomp,ggpubr,rstatix, ggtext)

       
setwd("~/Documents/ND_PhD/Writing Projects/Paper/MMAR_SecA2_ManuscriptData_2023/secA2_qPCR/")
data=data.frame(read.csv("secA2_3Biological_Replicates.csv",header=T))


#Make the strain levels a factor and order as desired
data$Strain <- factor(data$Strain, levels = c("WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"))

##Find the avaerage and sd of MARC3, RvC3, and delta_esxBA

mean_WT = mean(data$ΔΔCt[data$Strain=="WT"])
sd_WT = sd(data$ΔΔCt[data$Strain=="WT"])

mean_ΔsecA2 = mean(data$ΔΔCt[data$Strain=="ΔsecA2"])
sd_ΔsecA2 = sd(data$ΔΔCt[data$Strain=="ΔsecA2"])

mean_MARC3 = mean(data$ΔΔCt[data$Strain=="ΔsecA2/secA2_MM"])
sd_MARC3 = sd(data$ΔΔCt[data$Strain=="ΔsecA2/secA2_MM"])

mean_RvC3 = mean(data$ΔΔCt[data$Strain=="ΔsecA2/secA2_MT"])
sd_RvC3 = sd(data$ΔΔCt[data$Strain=="ΔsecA2/secA2_MT"])

mean_ΔesxBA = mean(data$ΔΔCt[data$Strain=="ΔesxBA"])
sd_ΔesxBA = sd(data$ΔΔCt[data$Strain=="ΔesxBA"])


########Bar plot with data points as dots#######
#colors
colors <- c("#990000", "#FFCC33", "#CC9900", "#996600", "#6699FF")

#Plot

L <- ggplot(data, aes(x=Strain, y=ΔΔCt, fill=Strain)) + geom_point(size=2) + theme_classic() + theme(axis.text=element_text(size=40),
                                                                                   axis.title.y=element_text(size=40),
                                                                                   axis.title.x=element_blank(),
                                                                                   axis.text.x=element_text(angle=45, hjust = 1)) 
#Add dots and bar
AveL <- L + stat_summary(geom = "bar", fun = "mean",color=c('#990000','#FFCC33','#CC9900','#996600','#6699FF'), linewidth=2, fill=c("white", "white","white","white", "white"), alpha=c(0.2)) + ylab(expression(paste(italic("secA2")," relative to ",italic("sigA")))) +
  scale_fill_manual('Strain', values=alpha(c('#990000','#FFCC33','#CC9900','#996600','#6699FF'), 0.1)) +
  theme(legend.position = "none", axis.text=element_text(size=40), 
        axis.title.y=element_text(size=40), axis.title.x=element_blank(), legend.title = element_blank(), legend.text = element_text(size=40), axis.text.x = element_markdown(angle=45, vjust=1, hjust=1),
        plot.title = element_text(hjust = 0, size=40)) +
  ###MUST USE ggtext with axis.text.x = element_markdown() and then any text to be in italics surrounded by * (ex. *secA2*) to get words/numbers italics in the text x-aixs labels (use the y version for y)
  scale_x_discrete(labels = c("WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*")) +
  ylab("Colony forming units \n (CFU)/mL at 2hpi")+ #Change the y-axis label
  coord_cartesian(ylim = c(0, 8))
 

#Manually add not-detected
PlotFinal <- AveL + geom_signif(annotation=c("nd", "***", "***", ""),
                       y_position = c(2,7.5,7.5,2),
                       xmin=c(2,3,4,5),
                       xmax=c(2,3,4,5), textsize = 14,
                       size=0) 

PlotFinal
ggsave("./secA2_CombinedData_BarPlot.png", width = 10, height = 10)




####STATISTICS????
###Use a Shapiro-Wilk normality test to determine if the data is normally distributed.
shapiro.test(data$ΔΔCt)
shapiro.test(data$ΔΔCt[data$Strain=="WT"])
shapiro.test(data$ΔΔCt[data$Strain=="ΔsecA2"])
shapiro.test(data$ΔΔCt[data$Strain=="ΔsecA2/secA2_MM"])
shapiro.test(data$ΔΔCt[data$Strain=="ΔsecA2/secA2_MT"])
shapiro.test(data$ΔΔCt[data$Strain=="ΔesxBA"])


##ALL Seem to be NORMAL (although the d_secA2 is very close!) so an ANVOA is okay?

##### ALTERNATIVELY TRY A KRUSKAL-WALLIS test
kruskal.test(ΔΔCt ~ Strain, data = data)

##As the p-value is less than 0.05, there is likely a significant interaction in the data set.
##Use a pair-wise comparison to determine which pair is significantly different.
##If comparing to WT, this shows the d_secA2 strain is statistically significant from WT. The others are not, although d_esxBA is close.
pairwise.wilcox.test(data$ΔΔCt, data$Strain, exact=F)


##Add One-Way ANOVA statistics
#compute One-way ANVOA
res.aov <- aov( ΔΔCt ~ Strain, data = data)

#Look at summary statistics
summary(res.aov)
#Can see there is statistical significance but can't see between which variables


#Run a Dunnett Test to determine which treatment are significant compared to the control (WT)
post_test <- glht(model = res.aov, linfct = mcp(Strain = "Dunnett"))


summary(post_test)
















#Run a Dunnett Test to determine which treatment are significant compared to the control (WT)
#Method 1

p_val <- DunnettTest(x=data$ΔΔCt, g=data$Strain) 
df <- p_val[["WT"]]
df <- df %>%
  as_tibble() %>%
  setNames(c('diff', 'lwr.ci', 'upr.ci', 'p_val'))

#Method2
post_test <- glht(res.aov,
                  linfct = mcp(Strain = "Dunnett")
)



df_p_val <- rstatix::t_test(data, ΔΔCt ~ Strain, ref.group = "WT") %>% 
  rstatix::add_xy_position()

compare_means(ΔΔCt ~ Strain,  data = data, ref.group = "WT",
              method = "anova")



#####Dot plot ONLY SD bars########



#Bar plot with no individual data points but standard deviation bars
x1 <- factor(D$Uninfected,levels=c("Uninfected", "WT", "delta_secA2", "delta_secA2/secA_2MM", "delta_secA2/secA2_MT", "delta_esxBA"))
D$strain <- x1

##Data and calculate SD, SE,and IC
DfDnew <- data %>% select(Strain, ΔΔCt)
mysum <- DfDnew %>%
  group_by(Strain) %>% 
  summarise(
    n=n(),
    mean= mean(ΔΔCt), 
    sd = sd(ΔΔCt)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))

#Plot

L <- ggplot(mysum, aes(x=Strain,y=mean)) + geom_point(size=0.5) + theme_classic() + theme(axis.text=element_text(size=5),
                                                                                          axis.title.y=element_text(size=5),
                                                                                          axis.title.x=element_blank(),
                                                                                          axis.text.x=element_text(angle=45, hjust = 1))
#Add dots and bar
AveL <- L + stat_summary(geom = "bar", fun = "mean", fill=c("#990066", rep("Black",5)), alpha=c(0.2)) + ylab("Relative Luminescence Unit (RLU)") +
  geom_errorbar( aes(x=Strain, ymin=mean-sd, ymax=mean+sd), width=0.2, colour="black", alpha=0.9, size=0.5)

AveL 

##Save plot
ggsave("/Users/lou/Documents/ND_PhD/Laboratory/Thesis/Data/Luminescence/L929-ISRE_PlateRead/2021.06.25/barplot_SDbars.png", width = 2, height = 2)

