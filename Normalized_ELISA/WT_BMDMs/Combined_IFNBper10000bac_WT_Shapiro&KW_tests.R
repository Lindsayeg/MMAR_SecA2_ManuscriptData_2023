
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(ggplot2,tidyverse,plyr,dplyr,DescTools,multcomp,ggpubr,rstatix,ggprism,magrittr,patchwork,ggtext,tibble)


#####ELISA
setwd("~/Documents/ND_PhD/Writing Projects/Paper/MMAR_SecA2_ManuscriptData_2023/Normalized_ELISA/WT_BMDMs/")
data=data.frame(read.csv("CD_WT_IFNB.csv",header=T, stringsAsFactors = FALSE))


#Removed uninfected and WT-HK
df <- data[c(1:18, 37:78, 97:138, 157:198, 217:258, 277:309, 316:336, 343:363, 370:390, 397:417, 424:432),c(1,2)]

#Replace ΔsecA2_MOI1 with ΔsecA2
df["Strain"][df["Strain"]=="ΔsecA2_MOI1"] <- "ΔsecA2"

#Make the strain levels a factor and order as desired
df$Strain <- factor(df$Strain, levels = c("Uninfected", "WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA", "WT-HK_MOI10"))



########Bar plot with data points as dots#######

#Plot IFNB concentration by strain

L <- ggplot(df, aes(x=Strain, y=IFNB_Concentration, fill=Strain)) + geom_point(size=0.1) + theme_classic() + theme(axis.text=element_text(size=16),
                                                                                                                   axis.title.y=element_text(size=16),
                                                                                                                   axis.title.x=element_blank(),
                                                                                                                   axis.text.x=element_text(angle=45, hjust = 1))
#Add data points as dots and bar graph for average values

AveL <- L + stat_summary(geom = "bar", fun = "mean", color=c("#999999", "#990000", "#FFCC33", "#CC9900","#996600", "#6699FF", "#999999"), size=0.4, fill=c("#999999", "#990000", "#FFCC33","#CC9900","#996600", "#6699FF", "#999999"), alpha=c(0.4)) + ylab("IFN-β Concentration \n (pg/ml)") +
  theme(axis.text.x=element_markdown(angle=45, hjust = 1, size=16), axis.title.x=element_blank(), legend.title = element_blank(), legend.position="none", 
        axis.title.y=element_text(size=16), plot.title = element_text(hjust = 0, size=16)) +
  #labs(title="IFN-β ELISA") +
  scale_x_discrete(labels = c("Uninfected"="Uninfected", "WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*", "WT-HK_MOI10" = "WT-HK_MOI10")) +
  coord_cartesian(ylim = c(0, 100))

AveL

###Use a Shapiro-Wilk normality test to determine if the data is normally distributed.

shapiro.test(df$IFNB_Concentration[df$Strain=="Uninfected"])
shapiro.test(df$IFNB_Concentration[df$Strain=="WT"])
shapiro.test(df$IFNB_Concentration[df$Strain=="ΔsecA2"])
shapiro.test(df$IFNB_Concentration[df$Strain=="ΔsecA2/secA2_MM"])
shapiro.test(df$IFNB_Concentration[df$Strain=="ΔsecA2/secA2_MT"])
shapiro.test(df$IFNB_Concentration[df$Strain=="ΔesxBA"])

###If it is normally distributed you can use a parametric test such as a one-way ANOVA
###to test for statistical significance. If it is not, use a non-parametric test 
###such as a Kruskal-Wallis test to examine statistical significance.
###You can also view the data by eye using hist(df$IFNB_Concentration), but the 
###shapiro test is more reliable.

##### ALTERNATIVELY TRY A KRUSKAL-WALLIS test
kruskal.test(IFNB_Concentration ~ Strain, data = df)

##As the p-value is less than 0.05, there is likely a significant interaction in the data set.
##Use a pair-wise comparison to determine which pair is significantly different.
##If comparing to WT, this shows the d_secA2 strain is statistically significant from WT. The others are not, although d_esxBA is close.
pairwise.wilcox.test(df$IFNB_Concentration, df$Strain, exact=F)

#Add statistics to graph
p1 <- AveL + geom_signif(annotation=c("***","ns","ns","ns","***","***"),
                         y_position = c(95,95,95,95,95,95),
                         xmin=c(1,3,4,5,6,7),
                         xmax=c(1,3,4,5,6,7), textsize = 6,
                         size=0) 

p1



###Find the average and sd
data_summary <- function(df, varname, groupnames){
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

Data_summary_ELISA <- data_summary(df, varname = "IFNB_Concentration", groupnames = c("Strain"))





######CFU
#Plot CFU by strain
setwd("~/Documents/ND_PhD/Laboratory/Thesis/Data/Excel/ELISA/CombinedData/Dose-response ELISA/WT BMDMS/")
data=data.frame(read.csv("CD_CFU.csv",header=T))

#Removed uninfected and WT-HK
df <- data[c(1:9, 19:39, 49:78, 97:138, 157:198, 217:249, 256:276, 283:303, 310:330, 337:357, 364:375),c(1,2)]

#Replace ΔsecA2_MOI1 with ΔsecA2
df["Strain"][df["Strain"]=="ΔsecA2_MOI1"] <- "ΔsecA2"
df["Strain"][df["Strain"]=="WT_HK"] <- "WT-HK_MOI10"


#Make the strain levels a factor and order as desired
df$Strain <- factor(df$Strain, levels = c("Uninfected", "WT", "ΔsecA2","ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA", "WT-HK_MOI10"))


L <- ggplot(df, aes(x=Strain, y=Uptake, fill=Strain)) + geom_point(size=0.1) + theme_classic() + theme(axis.text=element_text(size=16),
                                                                                                       axis.title.y=element_text(size=16),
                                                                                                       axis.title.x=element_blank(),
                                                                                                       axis.text.x=element_text(angle=45, hjust = 1))
#Add data points as dots and bar graph for average values
AveL <- L + stat_summary(geom = "bar", fun = "mean", color=c("#999999", "#990000", "#FFCC33","#CC9900","#996600", "#6699FF", "#999999"), linewidth=0.4, fill=c("#999999", "#990000", "#FFCC33","#CC9900","#996600", "#6699FF", "#999999"), alpha=c(0.8)) + ylab("Colony Forming Units \n (CFU/mL)") +
  theme(axis.text.x=element_markdown(angle=45, hjust = 1, size=16), axis.title.x=element_blank(), legend.title = element_blank(), legend.position="none",
        axis.title.y=element_text(size=16), plot.title = element_text(hjust = 0, size=16)) +
  scale_x_discrete(labels = c("Uninfected"="Uninfected", "WT" = "WT","ΔsecA2" = "Δ*secA2*","ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*", "WT-HK_MOI10" = "WT-HK_MOI10"))+
  #labs(title="No. Intracellular Bacteria at 2hpi") +
  coord_cartesian(ylim = c(0, 1000000))


AveL

###Use a Shapiro-Wilk normality test to determine if the data is normally distributed.
shapiro.test(df$Uptake)
shapiro.test(df$Uptake[df$Strain=="Uninfected"])
shapiro.test(df$Uptake[df$Strain=="WT"])
shapiro.test(df$Uptake[df$Strain=="ΔsecA2"])
shapiro.test(df$Uptake[df$Strain=="ΔsecA2/secA2_MM"])
shapiro.test(df$Uptake[df$Strain=="ΔsecA2/secA2_MT"])
shapiro.test(df$Uptake[df$Strain=="ΔesxBA"])

###If it is normally distributed you can use a parametric test such as a one-way ANOVA
###to test for statistical significance. If it is not, use a non-parametric test 
###such as a Kruskal-Wallis test to examine statistical significance.

##### ALTERNATIVELY TRY A KRUSKAL-WALLIS test
kruskal.test(Uptake ~ Strain, data = df)

##As the p-value is less than 0.05, there is likely a significant interaction in the data set.
##Use a pair-wise comparison to determine which pair is significantly different.
##If comparing to WT, this shows the d_secA2 strain is statistically significant from WT. The others are not, although d_esxBA is close.
pairwise.wilcox.test(df$Uptake, df$Strain, exact=F)


#Add statistics to graph
p2 <- AveL + geom_signif(annotation=c("***","***","***","***","***","***"),
                         y_position = c(900000,900000,900000,900000,900000,900000),
                         xmin=c(1,3,4,5,6,7),
                         xmax=c(1,3,4,5,6,7), textsize = 6,
                         size=0) 

p2


###Find the average and sd
data_summary <- function(df, varname, groupnames){
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

Data_summary_CFU <- data_summary(df, varname = "Uptake", groupnames = c("Strain"))

WT <- mean(df$Uptake[df$Strain=="WT"])
A <- mean(df$Uptake[df$Strain=="WT"])/mean(df$Uptake[df$Strain=="ΔsecA2"])
B <- mean(df$Uptake[df$Strain=="WT"])/mean(df$Uptake[df$Strain=="ΔsecA2/secA2_MM"])
C <- mean(df$Uptake[df$Strain=="WT"])/mean(df$Uptake[df$Strain=="ΔsecA2/secA2_MT"])
D <- mean(df$Uptake[df$Strain=="WT"])/mean(df$Uptake[df$Strain=="ΔesxBA"])

sum(A+B+C+D)/4



#######IFNBper10000bac.

setwd("~/Documents/ND_PhD/Laboratory/Thesis/Data/Excel/ELISA/CombinedData/Dose-response ELISA/WT BMDMS/")
data=data.frame(read.csv("Trial2_WTB.csv",header=T))


#Removed uninfected and WT-HK
df <- data[c(7:18, 37:54, 67:78, 97:114, 127:138, 157:174, 187:198, 217:234, 247:258, 277:294, 304:309, 316:324, 331:336, 343:351, 358:363, 370:378, 385:390, 397:405, 412:417, 424:432),c(1,2)]

#Replace ΔsecA2_MOI1 with ΔsecA2
df["Strain"][df["Strain"]=="ΔsecA2_MOI1"] <- "ΔsecA2"

#Make the strain levels a factor and order as desired
df$Strain <- factor(df$Strain, levels = c("WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"))



########Bar plot with data points as dots#######

#Plot IFNB concentration by strain

L <- ggplot(df, aes(x=Strain, y=Corrected, fill=Strain)) + geom_point(size=0.1) + theme_classic() + theme(axis.text=element_text(size=16),
                                                                                   axis.title.y=element_text(size=16),
                                                                                   axis.title.x=element_blank(),
                                                                                   axis.text.x=element_text(angle=45, hjust = 1))
#Add data points as dots and bar graph for average values
#AveL <- L + stat_summary(geom = "bar", fun = "mean", fill=c(rep("Black",5),"#990066"), alpha=c(0.2)) + ylab("IFN-β Concentration (pg/ml)") 
AveL <- L + stat_summary(geom = "bar", fun = "mean", color=c("#990000", "#FFCC33", "#CC9900","#996600", "#6699FF"), size=0.75, fill=c("#990000", "#FFCC33","#CC9900","#996600", "#6699FF"), alpha=c(0)) + ylab("IFN-β (pg/ml) \n per 10000 bac.") +
  theme(axis.text.x=element_markdown(angle=45, hjust = 1, size=16), axis.title.x=element_blank(), legend.title = element_blank(), legend.position="none", 
        axis.title.y=element_text(size=16), plot.title = element_text(hjust = 0, size=16)) +
  scale_x_discrete(labels = c("WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*")) +
  coord_cartesian(ylim = c(0, 15))

AveL

###Use a Shapiro-Wilk normality test to determine if the data is normally distributed.
shapiro.test(df$Corrected)
shapiro.test(df$Corrected[df$Strain=="WT"])
shapiro.test(df$Corrected[df$Strain=="ΔsecA2"])
shapiro.test(df$Corrected[df$Strain=="ΔsecA2/secA2_MM"])
shapiro.test(df$Corrected[df$Strain=="ΔsecA2/secA2_MT"])
shapiro.test(df$Corrected[df$Strain=="ΔesxBA"])

###If it is normally distributed you can use a parametric test such as a one-way ANOVA
###to test for statistical significance. If it is not, use a non-parametric test 
###such as a Kruskal-Wallis test to examine statistical significance.

#####KRUSKAL-WALLIS test
kruskal.test(Corrected ~ Strain, data = df)

##As the p-value is less than 0.05, there is likely a significant interaction in the data set.
##Use a pair-wise comparison to determine which pair is significantly different.
##If comparing to WT, this shows the d_secA2 strain is statistically significant from WT. The others are not, although d_esxBA is close.
pairwise.wilcox.test(df$Corrected, df$Strain, exact=F)

#Add statistics to graph
p3 <- AveL + geom_signif(annotation=c("***","***","***","***"),
                         y_position = c(14,14,14,14),
                         xmin=c(2,3,4,5),
                         xmax=c(2,3,4,5), textsize = 6,
                         size=0) 

p3

###Find the average and sd
data_summary <- function(df, varname, groupnames){
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

Data_summary_NormELISA <- data_summary(df, varname = "Corrected", groupnames = c("Strain"))






####Combine all figures into one
figure <- ggarrange(p2, p1, p3,
                    labels = c("A. ", "B. ", "C. "),
                    vjust = 27,
                    ncol = 3, nrow = 1)
figure

ggsave("./Combine_IFNBper10000bac.png", figure, width = 15, height = 5)







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

