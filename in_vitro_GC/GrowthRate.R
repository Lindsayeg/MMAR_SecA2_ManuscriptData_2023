###Data read into this script is an edited output from the CombinedData_GC_GrowthCurver script with the following modifications:
###(1) I changed the strains names to exclude the _1A, _1B, _1C endings.
###(2) I changed the growth rate values so they were separated by a period (.) and not a comma (,).
###This was done manually rather than through code.


## Install and load packages
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(ggplot2,tidyverse,plyr,dplyr,multcomp,ggpubr,reshape2,growthcurver,purrr,deSolve, ggtext,DescTools)

## Set working directory
setwd("~/Documents/ND_PhD/Laboratory/Thesis/Data/Excel/InVitro-GC/CombinedData/")

## Load data into tibble (tidyverse version of data frame - note that read_csv is the tidyverse version of read.csv
#file_name <- "~/Documents/ND_PhD/Laboratory/Thesis/Data/Excel/InVitro-GC/CombinedData/growthcurver_test.csv"
data = read.csv("GC_output_edited.csv", sep="\t", header=TRUE)


#Subset data
df <- data[,c("sample","r")]

#Make the strain levels a factor and order as desired
df$sample <- factor(df$sample, levels = c("WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"))


#Conduct an ANOVA
aov.results <- aov(r ~ sample, data=df)

#Look at summary statistics
#Can see there is statistical significance but can't see between which variables
summary(aov.results)


#Run a Dunnett Test to determine which treatment are significant compared to the control (WT)
post_test <- glht(model = aov.results, linfct = mcp(sample = "Dunnett"))

#Look at summary statistics
summary(post_test)


#Plot it

L <- ggplot(df, aes(x=sample, y=r)) + geom_point(size=0.1) + theme_classic() + theme(axis.text=element_text(size=12),
                                                                                                          axis.title.y=element_text(size=12),
                                                                                                          axis.title.x=element_blank(),
                                                                                                          axis.text.x=element_text(angle=45, hjust = 1, size=12))
#Add data points as dots and bar graph for average values
#AveL <- L + stat_summary(geom = "bar", fun = "mean", fill=c(rep("Black",5),"#990066"), alpha=c(0.2)) + ylab("IFN-β Concentration (pg/ml)") 
AveL <- L + stat_summary(geom = "bar", fun = "mean", color=c("#990000", "#FFCC33", "#CC9900","#996600", "#6699FF"), linewidth=0.4, fill=c("#990000", "#FFCC33","#CC9900","#996600", "#6699FF"), alpha=c(0)) + ylab("Bacterial growth rate (r)") + 
  theme(axis.text.x=element_markdown(angle=45, hjust = 1, size=12), axis.title.x=element_blank(), legend.title = element_blank(), legend.position="none", 
        axis.title.y=element_text(size=12), plot.title = element_text(hjust = 0, size=12)) +
  #labs(title="Growth rate (r)") +
  scale_x_discrete(labels = c("WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*", "ΔesxBA" = "Δ*esxBA*")) +
  coord_cartesian(ylim = c(0, 0.075))

AveL

#Manually add ANOVA results
p2 <- AveL + geom_signif(annotation=c("***", "", "", ""),
                         y_position = c(0.065, 0.065, 0.065, 0.065),
                         xmin=c(2,3,4,5),
                         xmax=c(2,3,4,5), textsize = 5,
                         size=0) 

p2


ggsave("./GrowrthRate_r.png", width = 4, height = 4)


###Use a Shapiro-Wilk normality test to determine if the data is normally distributed.
shapiro.test(df$r)
shapiro.test(df$r[df$sample=="WT"])
shapiro.test(df$r[df$sample=="ΔsecA2"])
shapiro.test(df$r[df$sample=="ΔsecA2/secA2_MM"])
shapiro.test(df$r[df$sample=="ΔsecA2/secA2_MT"])
shapiro.test(df$r[df$sample=="ΔesxBA"])

##### ALTERNATIVELY TRY A KRUSKAL-WALLIS test
kruskal.test(r ~ sample, data = df)

##As the p-value is less than 0.05, there is likely a significant interaction in the data set.
##Use a pair-wise comparison to determine which pair is significantly different.
##If comparing to WT, this shows the d_secA2 strain is statistically significant from WT. The others are not, although d_esxBA is close.
pairwise.wilcox.test(df$r, df$sample)







###DID NOT END UP USING THE BELOW CODE
##Data and calculate SD, SE,and IC
##Data and calculate SD, SE,and IC
#DfDnew <- df %>% select(sample, r)
mysum <- df %>%
  group_by(sample) %>% 
  dplyr::summarise(
    n=n(),
    mean= mean(r), 
    sd = sd(r)
  ) %>%
  mutate( se=sd/sqrt(n))  %>%
  mutate( ic=se * qt((1-0.05)/2 + .5, n-1))






