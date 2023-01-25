#########LINE GRAPH PLOT##########
## Install and load packages
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(ggplot2,tidyverse,plyr,dplyr,ggtext)

## Set working directory
setwd("~/Documents/ND_PhD/Writing Projects/Paper/MMAR_SecA2_ManuscriptData_2023/intracellular_GC/")

## Load data into tibble (tidyverse version of data frame - note that read_csv is the tidyverse version of read.csv
data = read_csv("CombinedData_COUNTS.csv")

#Remove NAs
data <- na.omit(data)


#colors
colors <- c("#990000", "#FFCC33", "#CC9900", "#996600", "#6699FF")


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
df2 <- data_summary(data, varname = "Count", groupnames = c("Strain", "Time"))
#Make the strain levels a factor and order as desired
df2$Strain <- factor(df2$Strain, levels = c("WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"))


## Plot it
P1 <- ggplot(df2, aes(x=Time, y=Count, color = Strain)) + theme_classic() +
  geom_point() +
  geom_line() +
  scale_color_manual(values=colors, labels = c("WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*")) +
  theme(legend.position = c(0.3,0.85), axis.text=element_text(size=15),
        axis.title.y=element_text(size=15), axis.title.x=element_text(size=15), legend.title = element_blank(), legend.text = element_markdown(size=15)) +
  ###MUST USE ggtext with axis.text.x = element_markdown() and then any text to be in italics surrounded by * (ex. *secA2*) to get words/numbers italics in the text x-aixs labels (use the y version for y)
  scale_y_log10() +
  ylab("Colony Forming Unit (CFU)/mL") + #Change the y-axis label +
  xlab("Hours post infection (hpi)") +   #Change the x-axis label
  geom_errorbar(aes(ymin=Count-sd, ymax=Count+sd), width=5, linewidth=0.5,
                position=position_dodge(1), color="black") 




P1

#ggsave("./LineGraph_LogScale.png", P1, width = 5, height = 5)


#######Uptake Graph

## Set working directory
setwd("~/Documents/ND_PhD/Writing Projects/Paper/MMAR_SecA2_ManuscriptData_2023/intracellular_GC/")

## Load data into tibble (tidyverse version of data frame - note that read_csv is the tidyverse version of read.csv
data = read_csv("CombinedData_Uptake.csv")

#Remove NAs
data <- na.omit(data)

##Find the fold change in bacterial counts of WT versus d_secA2 for bacterial uptake at 2hpi
mean(data$Count[data$Strain=="WT"])/mean(data$Count[data$Strain=="ΔsecA2"])


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
df3 <- data_summary(data, varname = "Count", groupnames = c("Strain"))
#Make the strain levels a factor and order as desired
df3$Strain <- factor(df3$Strain, levels = c("WT", "ΔsecA2", "ΔsecA2/secA2_MM", "ΔsecA2/secA2_MT", "ΔesxBA"))




## Plot it
P2 <- ggplot(df3, aes(x=Strain, y=Count, fill=Strain)) + theme_classic() +
  stat_summary(geom = "bar", fun = mean, position = "dodge", color=alpha(c('#990000','#FFCC33','#CC9900','#996600','#6699FF'), 1), linewidth=0.7) +
  scale_fill_manual('Strain', values=alpha(c('#990000','#FFCC33','#CC9900','#996600','#6699FF'), 0.1)) +
  theme(legend.position = "none", axis.text=element_text(size=15), 
        axis.title.y=element_text(size=15), axis.title.x=element_blank(), legend.title = element_blank(), legend.text = element_text(size=15), axis.text.x = element_markdown(angle=45, vjust=1, hjust=1),
        plot.title = element_text(hjust = 0, size=15)) +
  ###MUST USE ggtext with axis.text.x = element_markdown() and then any text to be in italics surrounded by * (ex. *secA2*) to get words/numbers italics in the text x-aixs labels (use the y version for y)
  scale_x_discrete(labels = c("WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/secA2_MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/secA2_MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*")) +
  ylab("Colony forming units \n (CFU)/mL at 2hpi")+ #Change the y-axis label
  #labs(title="Uptake") +
  geom_errorbar(aes(ymin=Count-sd, ymax=Count+sd), width=0.2, size=0.5, color="black")+
  coord_cartesian(ylim = c(0, 40000))


###Use a Shapiro-Wilk normality test to determine if the data is normally distributed.
shapiro.test(data$Count)
shapiro.test(data$Count[data$Strain=="WT"])
shapiro.test(data$Count[data$Strain=="ΔsecA2"])
shapiro.test(data$Count[data$Strain=="ΔsecA2/secA2_MM"])
shapiro.test(data$Count[data$Strain=="ΔsecA2/secA2_MT"])
shapiro.test(data$Count[data$Strain=="ΔesxBA"])


##### ALTERNATIVELY TRY A KRUSKAL-WALLIS test
kruskal.test(Count ~ Strain, data = data)

##As the p-value is less than 0.05, there is likely a significant interaction in the data set.
##Use a pair-wise comparison to determine which pair is significantly different.
##If comparing to WT, this shows the d_secA2 strain is statistically significant from WT. The others are not, although d_esxBA is close.
pairwise.wilcox.test(data$Count, data$Strain, exact=F)

#Add statistics to graph
P2 <- P2 + geom_signif(annotation=c("***","***","***","***"),
                         y_position = c(34000,34000,34000,34000),
                         xmin=c(2,3,4,5),
                         xmax=c(2,3,4,5), textsize = 6,
                         size=0) 

P2

## Save plot
#ggsave("./2022.07.06-10_Uptake.png", P, width = 5, height = 5)



####
###
##
#Combine all figures into one
figure <- ggarrange(P1, P2,
                    labels = c("A. ", "B. "),
                    vjust = 33,
                    ncol = 2, nrow = 1)
figure

ggsave("./Line&Uptake_Graphs.png", width = 10, height = 5)





