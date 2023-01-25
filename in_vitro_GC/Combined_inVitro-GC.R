## Install and load packages
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(ggplot2,tidyverse,plyr,dplyr,ggtext)

## Set working directory
setwd("~/Documents/ND_PhD/Writing Projects/Paper/MMAR_SecA2_ManuscriptData_2023/in_vitro_GC/")

## Load data into tibble (tidyverse version of data frame - note that read_csv is the tidyverse version of read.csv
data = data.frame(read_csv("Combined-inVitro-GC.csv"))

#Replace ΔsecA2(8) with ΔsecA2
data["Strain"][data["Strain"]=="ΔsecA2(8)"] <- "ΔsecA2"

#Change column names
#names(data) <- factor(c("Strain", "0", "24", "48", "72", "96", "120"))
names(data) <- factor(c("Strain", "0", "24", "48", "72", "96", "120", "144", "168"))

##Make the data into a long format
data_long <- gather(data, Time, OD, as.character(seq(0,168,24)), factor_key=TRUE)
#data_long <- gather(data, Time, OD, as.character(seq(0,120,24)), factor_key=TRUE)


##Make the time data numeric
data_long$Time <- as.numeric(as.character(data_long$Time))


#Make the strain levels a factor and order as desired
data_long$Strain <- factor(data_long$Strain, levels = c("WT", "ΔsecA2", "ΔsecA2/psecA2MM", "ΔsecA2/psecA2MT", "ΔesxBA"))

#Remove NA values
data_long <- na.omit(data_long)

#colors
colors <- c("#990000", "#FFCC33","#CC9900","#996600", "#6699FF")

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
df2 <- data_summary(data_long, varname = "OD", groupnames = c("Strain", "Time"))
#Make the strain levels a factor and order as desired
df2$Strain <- factor(df2$Strain, levels = c("WT", "ΔsecA2", "ΔsecA2/psecA2MM", "ΔsecA2/psecA2MT", "ΔesxBA"))


## Plot it
P <- ggplot(df2, aes(x=Time, y=OD, color = Strain)) + theme_classic() +
  geom_point() +
  geom_line() +
  scale_color_manual(values=colors, labels = c("WT" = "WT","ΔsecA2" = "Δ*secA2*", "ΔsecA2/psecA2MM" = "Δ*secA2*/p*secA2<sub><i>MM</i></sub>*","ΔsecA2/psecA2MT" = "Δ*secA2*/p*secA2<sub><i>MT</i></sub>*", "ΔesxBA" = "Δ*esxBA*")) +
  theme(legend.position = c(0.72,0.25), axis.text=element_text(size=12),
        axis.title.y=element_text(size=12), axis.title.x=element_text(size=12), legend.title = element_blank(), legend.text = element_markdown(size=12)) +
  scale_y_log10() +
  ylab("Optical Density (OD600)") + #Change the y-axis label +
  xlab("Time (hours)") +   #Change the x-axis label
  geom_errorbar(aes(ymin=OD-sd, ymax=OD+sd), width=5, linewidth=0.5,
                position=position_dodge(1), color="black") 



P

## Save plot
ggsave("./Triplicate.png", P, width = 4, height = 4)






