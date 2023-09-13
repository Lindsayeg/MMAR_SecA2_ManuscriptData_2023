## Install and load packages
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(ggplot2,tidyverse,plyr,dplyr,ggtext,rstatix,ggpubr,DescTools)

## Set working directory
setwd("~/Desktop/")

## Load data into tibble (tidyverse version of data frame - note that read_csv is the tidyverse version of read.csv
data = read_csv("qPRC_8hpi_BMDMs.csv")

#Remove NAs
data <- na.omit(data)


#+++++++++++++++++++++++++
# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable
#to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
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

###
df2 <- data_summary(data, varname = "Value", groupnames = c("Transcript", "Strain"))

#Make the strain levels a factor and order as desired
df2$Strain <- factor(df2$Strain, levels = c("Uninfected", "WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"))




gg <- ggplot(df2, aes(x=Strain, y=Value, fill=Strain)) +
  theme_classic() +
  geom_bar(stat="identity", color="black", position=position_dodge()) +
  facet_wrap(~Transcript, scales = "free_y")

P <- gg + geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), size =0.2, width=0.4, position=position_dodge()) + facet_wrap(~Transcript, scales = "free_y") +
  scale_fill_manual('Strain',values=alpha(c('#990000','#FFCC33','#CC9900','#996600','#6699FF'), 0.9)) +
  theme(legend.position = "none", axis.text=element_text(size=11),
                 axis.title.y=element_text(size=11), axis.title.x=element_blank(), axis.text.x = element_markdown(angle=45, vjust=1, hjust=1, size=11)) +
  scale_x_discrete(labels = c("WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*")) +
  scale_y_log10() +
  ylab("Transcript abundance \n relative to GAPDH")

#Linear scale rather than log10
B <- gg + geom_errorbar(aes(ymin=Value-sd, ymax=Value+sd), size =0.2, width=0.4, position=position_dodge()) + facet_wrap(~Transcript, scales = "free_y") +
  scale_fill_manual('Strain',values=alpha(c('#990000','#FFCC33','#CC9900','#996600','#6699FF'), 0.9)) +
  theme(legend.position = "none", axis.text=element_text(size=11),
        axis.title.y=element_text(size=11), axis.title.x=element_blank(), axis.text.x = element_markdown(angle=45, vjust=1, hjust=1, size=11)) +
  scale_x_discrete(labels = c("WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*")) +
  ylab("Transcript abundance \n relative to GAPDH")


P
B

#####Test data for normality 

###Using a histogram
# histogram
res_aov <- aov(Value ~ Transcript,
               data = data
)
hist(res_aov$residuals)

# QQ-plot
library(car)
qqPlot(res_aov$residuals,
       id = FALSE # id = FALSE to remove point identification
)

##using a Shapiro-Wilks Test. 
###If p<0.05, then data is not normal and a non-parametric test is required.

shapiro.test(data$Value)
shapiro.test(data$Value[data$Transcript=="Ifn-β"])
shapiro.test(data$Value[data$Transcript=="Irf7"])
shapiro.test(data$Value[data$Transcript=="Rig-I"])



####For non-parametric data, use a Kruskal-wallis test

kruskal.test(Value ~ Transcript, data = data)


####If p < 0.05 then there is a significant interaction. Do a pairwise comparison.
#StatsOutpout <- data %>%
# group_by(Strain)%>%
# wilcox_test(Geometric.Mean ~ Sera.Ab, ref.group = "NMS") %>%
# adjust_pvalue(method = "bonferroni") %>%
# add_significance("p.adj")



###Anova by trancript 
df_ifnb <- data[c(1:12, 31:40), c(1,2)]
#Set Control - here I want everything compared to WT and not Uninfected so WT must be listed first
df_ifnb$Strain <- factor(df_ifnb$Strain, levels = c("WT","ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"), ordered = TRUE)
#compute One-way ANVOA
res.aov <- aov(Value ~ Strain, data = df_ifnb)
DunnettTest(Value ~ Strain, df_ifnb)

df_irf7 <- data[c(13:24, 41:50), c(1,2)]
#Set Control - here I want everything compared to WT and not Uninfected so WT must be listed first
df_irf7$Strain <- factor(df_irf7$Strain, levels = c("WT","ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"), ordered = TRUE)
#compute One-way ANVOA
res.aov <- aov(Value ~ Strain, data = df_irf7)
DunnettTest(Value ~ Strain, df_irf7)

df_rigI <- data[c(25:36, 51:60), c(1,2)]
#Set Control - here I want everything compared to WT and not Uninfected so WT must be listed first
df_rigI$Strain <- factor(df_rigI$Strain, levels = c("WT","ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"), ordered = TRUE)
#compute One-way ANVOA
res.aov <- aov(Value ~ Strain, data = df_rigI)
DunnettTest(Value ~ Strain, df_rigI)




###Are these the right statistics?
###Statistics by transcript group using rStatix

df_ifnb <- data[c(1:12, 31:40), c(1,2)]

stat.test_IFNb <- df_ifnb %>%
  t_test(Value ~ Strain, ref.group = "WT") 
stat.test_IFNb

df_irf7 <- data[c(13:24, 41:50), c(1,2)]

stat.test_Irf7 <- df_irf7 %>%
  t_test(Value ~ Strain, ref.group = "WT") 
stat.test_Irf7

df_rigI <- data[c(25:36, 51:60), c(1,2)]

stat.test_rigI <- df_rigI %>%
  t_test(Value ~ Strain, p.adjust.method="bonferroni", ref.group = "WT") 
stat.test_rigI



#Statistics using rStatix
stat.test <- data %>%
  group_by(Transcript) %>%
  t_test(Value ~ Strain, p.adjust.method="bonferroni", ref.group = "WT") 
stat.test



###Plot it with statistics. Here ** is p<0.01  
Q <- P + geom_signif(annotations = c('','','', '**'),
                                 y_position = c(1,1,1,1),
                                 xmin=c(1,1,1,1),
                                 xmax=c(2,3,4,5), textsize = 4,
                                 tip_length=0, vjust=0.5)

R <- B + geom_signif(annotations = c('','','', '**'),
                     y_position = c(5,5,5,5),
                     xmin=c(1,1,1,1),
                     xmax=c(2,3,4,5), textsize = 4,
                     tip_length=0, vjust=0.5)


Q
R
## Save plot
ggsave("./Faceted_data.png", Q, width = 5, height = 5)
ggsave("./Faceted_data_R.png", R, width = 5, height = 5)








