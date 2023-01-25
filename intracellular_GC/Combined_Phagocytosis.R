## Install and load packages
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(ggplot2,tidyverse,plyr,dplyr,ggtext)

## Set working directory
setwd("~/Documents/ND_PhD/Writing Projects/Paper/MMAR_SecA2_ManuscriptData_2023/intracellular_GC/")

## Load data into tibble (tidyverse version of data frame - note that read_csv is the tidyverse version of read.csv
data = read_csv("CombinedData_Phagocytosis.csv")

#Remove NAs
data <- na.omit(data)


#colors
colors <- c("#990000", "#FFCC33","#CC9900","#996600", "#6699FF","#990000", "#FFCC33","#CC9900","#996600", "#6699FF")
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
df2 <- data_summary(data, varname = "Rel.Uptake", groupnames = c("Strain"))
#Make the strain levels a factor and order as desired
df2$Strain <- factor(df2$Strain, levels = c("WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"))




## Plot it
P <- ggplot(df2, aes(x=Strain, y=Rel.Uptake, fill=Strain)) + theme_classic() +
  stat_summary(geom = "bar", fun = mean, position = "dodge", color=alpha(c('#990000','#FFCC33','#CC9900','#996600','#6699FF'), 1), linewidth=0.7) +
  scale_fill_manual('Strain', values=alpha(c('#990000','#FFCC33','#CC9900','#996600','#6699FF'), 0.6)) +
  theme(legend.position = "none", axis.text=element_text(size=8), 
        axis.title.y=element_text(size=8), axis.title.x=element_blank(), legend.title = element_blank(), legend.text = element_text(size=8), axis.text.x = element_markdown(angle=45, vjust=1, hjust=1),
        plot.title = element_text(hjust = 0, size=10)) +
  ###MUST USE ggtext with axis.text.x = element_markdown() and then any text to be in italics surrounded by * (ex. *secA2*) to get words/numbers italics in the text x-aixs labels (use the y version for y)
  scale_x_discrete(labels = c("WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*")) +
  ylab("Bacterial uptake relative to input")+ #Change the y-axis label
  labs(title="Input Relative to Uptake") +
  geom_errorbar(aes(ymin=Rel.Uptake-sd, ymax=Rel.Uptake+sd), width=0.2, size=0.5, color="black") +
  coord_cartesian(ylim = c(0, 1.5))





P



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
shapiro.test(data$Rel.Uptake)
shapiro.test(data$Rel.Uptake[data$Strain=="WT"])
shapiro.test(data$Rel.Uptake[data$Strain=="ΔsecA2"])
shapiro.test(data$Rel.Uptake[data$Strain=="ΔsecA2/secA2_MM"])
shapiro.test(data$Rel.Uptake[data$Strain=="ΔsecA2/secA2_MT"])
shapiro.test(data$Rel.Uptake[data$Strain=="ΔesxBA"])


##Add One-Way ANOVA statistics
#Set Control - here I want everything compared to WT and not Uninfected so WT must be listed first
data$Strain <- factor(data$Strain, levels = c("WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"), ordered = TRUE)
#compute One-way ANVOA
res.aov <- aov(Rel.Uptake ~ Strain, data = data)

#Look at summary statistics
#Can see there is statistical significance but can't see between which variables
summary(res.aov)


#Run a Dunnett Test to determine which treatment are significant compared to the control (WT)
post_test <- glht(model = res.aov, linfct = mcp(Strain = "Dunnett"))

#Look at summary statistics
summary(post_test)



#Manually add ANOVA results (SEE BELOW)
p2 <- P + geom_signif(annotation=c("ns", "ns", "ns", "ns"),
                         y_position = c(1.1,1.2,1.3,1.4),
                         xmin=c(2,3,4,5),
                         xmax=c(1,1,1,1), textsize = 2) 

p2


## Save plot
ggsave("./CD_Phagocytosis.png", p2, width = 2.5, height = 2.5)

