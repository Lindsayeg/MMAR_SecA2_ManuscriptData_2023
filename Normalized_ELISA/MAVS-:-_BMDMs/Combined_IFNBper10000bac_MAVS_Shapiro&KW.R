
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(ggplot2,tidyverse,plyr,dplyr,DescTools,multcomp,ggpubr,rstatix,ggprism,magrittr,patchwork,ggtext,tibble)


#####ELISA
setwd("~/Documents/ND_PhD/Writing Projects/Paper/MMAR_SecA2_ManuscriptData_2023/Normalized_ELISA/MAVS-:-_BMDMs/")
data=data.frame(read.csv("CD_MAVS_IFNB.csv",header=T, stringsAsFactors = FALSE))


#Make the strain levels a factor and order as desired
data$Strain <- factor(data$Strain, levels = c("Uninfected", "WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA", "WT-HK_MOI10"))



########Bar plot with data points as dots#######

#Plot IFNB concentration by strain

L <- ggplot(data, aes(x=Strain, y=IFNB_Concentration, fill=Strain)) + geom_point(size=0.1) + theme_classic() + theme(axis.text=element_text(size=14),
                                                                                                                   axis.title.y=element_text(size=14),
                                                                                                                   axis.title.x=element_blank(),
                                                                                                                   axis.text.x=element_text(angle=45, hjust = 1))
#Add data points as dots and bar graph for average values

AveL <- L + stat_summary(geom = "bar", fun = "mean", color=c("#999999", "#990000", "#FFCC33", "#CC9900","#996600", "#6699FF", "#999999"), size=0.4, fill=c("#999999", "#990000", "#FFCC33","#CC9900","#996600", "#6699FF", "#999999"), alpha=c(0.4)) + ylab("IFN-β Concentration \n (pg/ml)") +
  theme(axis.text.x=element_markdown(angle=45, hjust = 1, size=14), axis.title.x=element_blank(), legend.title = element_blank(), legend.position="none", 
        axis.title.y=element_text(size=14), plot.title = element_text(hjust = 0, size=16)) +
  #labs(title="IFN-β ELISA") +
  scale_x_discrete(labels = c("Uninfected"="Uninfected", "WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*", "WT-HK_MOI10" = "WT-HK_MOI10")) +
  coord_cartesian(ylim = c(0, 125))

AveL

###Use a Shapiro-Wilk normality test to determine if the data is normally distributed.
shapiro.test(data$IFNB_Concentration)
shapiro.test(data$IFNB_Concentration[data$Strain=="WT"])
shapiro.test(data$IFNB_Concentration[data$Strain=="ΔsecA2"])
shapiro.test(data$IFNB_Concentration[data$Strain=="ΔsecA2/secA2_MM"])
shapiro.test(data$IFNB_Concentration[data$Strain=="ΔsecA2/secA2_MT"])
shapiro.test(data$IFNB_Concentration[data$Strain=="ΔesxBA"])

###If it is normally distributed you can use a parametric test such as a one-way ANOVA
###to test for statistical significance. If it is not, use a non-parametric test 
###such as a Kruskal-Wallis test to examine statistical significance.
###You can also view the data by eye using hist(df$IFNB_Concentration), but the 
###shapiro test is more reliable.

##### ALTERNATIVELY TRY A KRUSKAL-WALLIS test
kruskal.test(IFNB_Concentration ~ Strain, data = data)

##As the p-value is less than 0.05, there is likely a significant interaction in the data set.
##Use a pair-wise comparison to determine which pair is significantly different.
##If comparing to WT, this shows the d_secA2 strain is statistically significant from WT. The others are not, although d_esxBA is close.
pairwise.wilcox.test(data$IFNB_Concentration, data$Strain, exact=F)

#Add statistics to graph
p1 <- AveL + geom_signif(annotation=c("***","ns","ns","ns","***","***"),
                         y_position = c(115,115,115,115,115,115),
                         xmin=c(1,3,4,5,6,7),
                         xmax=c(1,3,4,5,6,7), textsize = 6,
                         size=0) 

p1


###Find the average and sd
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

Data_summary_ELISA <- data_summary(data, varname = "IFNB_Concentration", groupnames = c("Strain"))







######CFU
#Plot CFU by strain
setwd("~/Documents/ND_PhD/Writing Projects/Paper/MMAR_SecA2_ManuscriptData_2023/Normalized_ELISA/MAVS-:-_BMDMs/")
data=data.frame(read.csv("CD_MAVS_CFU.csv",header=T))

#Removed uninfected and WT-HK
df <- data[c(1:9, 16:36, 43:63, 70:90, 97:117, 124:135), ]
#Remove the third column
data <- df[,c(1:2)]

#Replace ΔsecA2_MOI1 with ΔsecA2
data["Strain"][data["Strain"]=="ΔsecA2_MOI1"] <- "ΔsecA2"

#Make the strain levels a factor and order as desired
data$Strain <- factor(data$Strain, levels = c("Uninfected", "WT", "ΔsecA2","ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA", "WT-HK_MOI10"))


L <- ggplot(data, aes(x=Strain, y=Uptake, fill=Strain)) + geom_point(size=0.1) + theme_classic() + theme(axis.text=element_text(size=14),
                                                                                                       axis.title.y=element_text(size=14),
                                                                                                       axis.title.x=element_blank(),
                                                                                                       axis.text.x=element_text(angle=45, hjust = 1))
#Add data points as dots and bar graph for average values
AveL <- L + stat_summary(geom = "bar", fun = "mean", color=c("#999999", "#990000", "#FFCC33","#CC9900","#996600", "#6699FF", "#999999"), linewidth=0.4, fill=c("#999999", "#990000", "#FFCC33","#CC9900","#996600", "#6699FF", "#999999"), alpha=c(0.8)) + ylab("Colony Forming Units \n (CFU/mL)") +
  theme(axis.text.x=element_markdown(angle=45, hjust = 1, size=14), axis.title.x=element_blank(), legend.title = element_blank(), legend.position="none",
        axis.title.y=element_text(size=14), plot.title = element_text(hjust = 0, size=16)) +
  scale_x_discrete(labels = c("Uninfected"="Uninfected", "WT" = "WT","ΔsecA2" = "Δ*secA2*","ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*", "WT-HK_MOI10" = "WT-HK_MOI10"))+
  #labs(title="No. Intracellular Bacteria at 2hpi") +
  coord_cartesian(ylim = c(0, 600000))


AveL

###Use a Shapiro-Wilk normality test to determine if the data is normally distributed.
shapiro.test(data$Uptake)
shapiro.test(data$Uptake[data$Strain=="WT"])
shapiro.test(data$Uptake[data$Strain=="ΔsecA2"])
shapiro.test(data$Uptake[data$Strain=="ΔsecA2/secA2_MM"])
shapiro.test(data$Uptake[data$Strain=="ΔsecA2/secA2_MT"])
shapiro.test(data$Uptake[data$Strain=="ΔesxBA"])

###If it is normally distributed you can use a parametric test such as a one-way ANOVA
###to test for statistical significance. If it is not, use a non-parametric test 
###such as a Kruskal-Wallis test to examine statistical significance.

##### ALTERNATIVELY TRY A KRUSKAL-WALLIS test
kruskal.test(Uptake ~ Strain, data = data)

##As the p-value is less than 0.05, there is likely a significant interaction in the data set.
##Use a pair-wise comparison to determine which pair is significantly different.
##If comparing to WT, this shows the d_secA2 strain is statistically significant from WT. The others are not, although d_esxBA is close.
pairwise.wilcox.test(data$Uptake, data$Strain, exact=F)


#Add statistics to graph
p2 <- AveL + geom_signif(annotation=c("***","*","*","**","**","***"),
                         y_position = c(550000,550000,550000,550000,550000,550000),
                         xmin=c(1,3,4,5,6,7),
                         xmax=c(1,3,4,5,6,7), textsize = 6,
                         size=0) 

p2

#Make the strain levels a factor and order as desired
data$Strain <- factor(data$Strain, levels = c("WT", "Uninfected", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA", "WT-HK_MOI10"))

#Conduct an ANOVA
aov.results <- aov(Uptake ~ Strain, data=data)

#Look at summary statistics
#Can see there is statistical significance but can't see between which variables
summary(aov.results)


#Run a Dunnett Test to determine which treatment are significant compared to the control (WT)
post_test <- glht(model = aov.results, linfct = mcp(Strain = "Dunnett"))

#Look at summary statistics
summary(post_test)

###Find the average and sd
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

Data_summary_CFU <- data_summary(data, varname = "Uptake", groupnames = c("Strain"))

WT <- mean(data$Uptake[data$Strain=="WT"])
A <- mean(data$Uptake[data$Strain=="WT"])/mean(data$Uptake[data$Strain=="ΔsecA2"])
B <- mean(data$Uptake[data$Strain=="WT"])/mean(data$Uptake[data$Strain=="ΔsecA2/secA2_MM"])
C <- mean(data$Uptake[data$Strain=="WT"])/mean(data$Uptake[data$Strain=="ΔsecA2/secA2_MT"])
D <- mean(data$Uptake[data$Strain=="WT"])/mean(data$Uptake[data$Strain=="ΔesxBA"])

sum(A+B+C+D)/4


#######IFNBper10000bac.

setwd("~/Documents/ND_PhD/Writing Projects/Paper/MMAR_SecA2_ManuscriptData_2023/Normalized_ELISA/MAVS-:-_BMDMs/")
data=data.frame(read.csv("CD_MAVS_IFNBper10000Bac.csv",header=T))


#Removed uninfected and WT-HK
df <- data[c(4:9, 16:24, 31:36, 43:51, 58:63, 70:78, 85:90, 97:105, 112:117, 124:132),c(1,2)]

#Replace ΔsecA2_MOI1 with ΔsecA2
df["Strain"][df["Strain"]=="ΔsecA2_MOI1"] <- "ΔsecA2"

#Make the strain levels a factor and order as desired
df$Strain <- factor(df$Strain, levels = c("WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"))



########Bar plot with data points as dots#######

#Plot IFNB concentration by strain

L <- ggplot(df, aes(x=Strain, y=IFNBper10000bac, fill=Strain)) + geom_point(size=0.1) + theme_classic() + theme(axis.text=element_text(size=14),
                                                                                                                axis.title.y=element_text(size=14),
                                                                                   axis.title.x=element_blank(),
                                                                                   axis.text.x=element_text(angle=45, hjust = 1))
#Add data points as dots and bar graph for average values
#AveL <- L + stat_summary(geom = "bar", fun = "mean", fill=c(rep("Black",5),"#990066"), alpha=c(0.2)) + ylab("IFN-β Concentration (pg/ml)") 
AveL <- L + stat_summary(geom = "bar", fun = "mean", color=c("#990000", "#FFCC33", "#CC9900","#996600", "#6699FF"), size=0.75, fill=c("#990000", "#FFCC33","#CC9900","#996600", "#6699FF"), alpha=c(0)) + ylab("IFN-β (pg/ml) per \n 10000 bac.") +
  theme(axis.text.x=element_markdown(angle=45, hjust = 1, size=14), axis.title.x=element_blank(), legend.title = element_blank(), legend.position="none", 
        axis.title.y=element_text(size=14), plot.title = element_text(hjust = 0, size=16)) +
  scale_x_discrete(labels = c("WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*")) +
  coord_cartesian(ylim = c(0, 20))

AveL

###Use a Shapiro-Wilk normality test to determine if the data is normally distributed.
shapiro.test(df$IFNBper10000bac)
shapiro.test(df$IFNBper10000bac[df$Strain=="WT"])
shapiro.test(df$IFNBper10000bac[df$Strain=="ΔsecA2"])
shapiro.test(df$IFNBper10000bac[df$Strain=="ΔsecA2/secA2_MM"])
shapiro.test(df$IFNBper10000bac[df$Strain=="ΔsecA2/secA2_MT"])
shapiro.test(df$IFNBper10000bac[df$Strain=="ΔesxBA"])

###If it is normally distributed you can use a parametric test such as a one-way ANOVA
###to test for statistical significance. If it is not, use a non-parametric test 
###such as a Kruskal-Wallis test to examine statistical significance.

##### ALTERNATIVELY TRY A KRUSKAL-WALLIS test
kruskal.test(IFNBper10000bac ~ Strain, data = df)

##As the p-value is less than 0.05, there is likely a significant interaction in the data set.
##Use a pair-wise comparison to determine which pair is significantly different.
##If comparing to WT, this shows the d_secA2 strain is statistically significant from WT. The others are not, although d_esxBA is close.
pairwise.wilcox.test(df$IFNBper10000bac, df$Strain, exact=F)

#Add statistics to graph
p3 <- AveL + geom_signif(annotation=c("ns","**","**","*"),
                         y_position = c(19,19,19,19),
                         xmin=c(2,3,4,5),
                         xmax=c(2,3,4,5), textsize = 6,
                         size=0) 

p3





####Combine all figures into one
figure <- ggarrange(p2, p1, p3,
                    labels = c("A. ", "B. ", "C. "),
                    vjust = 27,
                    ncol = 3, nrow = 1)
figure

ggsave("./Combine_MAVS.png", width = 12, height = 4)







########ANOVA Example
####
##Add One-Way ANOVA statistics
#Set Control - here I want everything compared to WT and not Uninfected so WT must be listed first
data$Strain <- factor(data$Strain, levels = c("WT", "Uninfected", "ΔsecA2","ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA", "WT-HK_MOI10"), ordered = TRUE)
#compute One-way ANVOA
res.aov <- aov(Corrected ~ Strain, data = df)

#Look at summary statistics
#Can see there is statistical significance but can't see between which variables
summary(res.aov)


#Run a Dunnett Test to determine which treatment are significant compared to the control (WT)
post_test <- glht(model = res.aov, linfct = mcp(Strain = "Dunnett"))

#Look at summary statistics
summary(post_test)

