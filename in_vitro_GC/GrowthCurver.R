## Install and load packages
if(!require(pacman)){install.packages("pacman");library(pacman)}
p_load(ggplot2,tidyverse,plyr,dplyr,multcomp,ggpubr,reshape2,growthcurver,purrr,deSolve)

## Set working directory
setwd("~/Documents/ND_PhD/Writing Projects/Paper/MMAR_SecA2_ManuscriptData_2023/in_vitro_GC/")

## Load data into tibble (tidyverse version of data frame - note that read_csv is the tidyverse version of read.csv
#file_name <- "~/Documents/ND_PhD/Laboratory/Thesis/Data/Excel/InVitro-GC/CombinedData/growthcurver_test.csv"
data = read_csv("Growthcurver.csv")


#Open png file
png(filename="Growthcurvere_3BioReps_Equations.png", width=12, height=8, unit="in", pointsize = 12, res=900)


#Make the graph so that it combines 15 plots in a 5x3 matrix
par(mfrow=c(3,5),mar=c(2.1,2.1,3.1,0.1),oma=c(2,2,0,1))


### Fit data to curve using growthcurver

gc_fit_WT_A1 <- SummarizeGrowth(data$Time, data$A1)
gc_fit_ΔsecA2_B1 <- SummarizeGrowth(data$Time, data$B1)
gc_fit_MMAR_C1 <- SummarizeGrowth(data$Time, data$C1)
gc_fit_Rv_D1 <- SummarizeGrowth(data$Time, data$D1)
gc_fit_ΔesxBA_E1 <- SummarizeGrowth(data$Time, data$E1)

gc_fit_WT_A2 <- SummarizeGrowth(data$Time, data$A2)
gc_fit_ΔsecA2_B2 <- SummarizeGrowth(data$Time, data$B2)
gc_fit_MMAR_C2 <- SummarizeGrowth(data$Time, data$C2)
gc_fit_Rv_D2 <- SummarizeGrowth(data$Time, data$D2)
gc_fit_ΔesxBA_E2 <- SummarizeGrowth(data$Time, data$E2)

gc_fit_WT_A3 <- SummarizeGrowth(data$Time, data$A3)
gc_fit_ΔsecA2_B3 <- SummarizeGrowth(data$Time, data$B3)
gc_fit_MMAR_C3 <- SummarizeGrowth(data$Time, data$C3)
gc_fit_Rv_D3 <- SummarizeGrowth(data$Time, data$D3)
gc_fit_ΔesxBA_E3 <- SummarizeGrowth(data$Time, data$E3)


##Plot Data

plot(data$Time, data$A1, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(paste("WT_1: r=0.044",sep="")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_WT_A1$model),
           col="red",
           lwd=2))

plot(data$Time, data$B1, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(paste("Δ",italic("secA2"),"_1: r=0.032",sep="")), ylim=c(0,12), 
     lines(data$Time,
           predict(gc_fit_ΔsecA2_B1$model),
           col="red",
           lwd=2))

plot(data$Time, data$C1, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(atop(paste("Δ",italic("secA2"),"/p",italic("secA2"[MM]),"_1:",sep=""), "r=0.044")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_MMAR_C1$model),
           col="red",
           lwd=2))

plot(data$Time, data$D1, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(atop(paste("Δ",italic("secA2"),"/p",italic("secA2"[MT]),"_1:",sep=""), "r=0.043")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_Rv_D1$model),
           col="red",
           lwd=2))

plot(data$Time, data$E1, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(paste("Δ",italic("esxBA"),"_1: r=0.041",sep="")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_ΔesxBA_E1$model),
           col="red",
           lwd=2))

plot(data$Time, data$A2, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(paste("WT_2: r=0.043",sep="")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_WT_A2$model),
           col="red",
           lwd=2))

plot(data$Time, data$B2, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(paste("Δ",italic("secA2"),"_2: r=0.029",sep="")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_ΔsecA2_B2$model),
           col="red",
           lwd=2))

plot(data$Time, data$C2, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(atop(paste("Δ",italic("secA2"),"/p",italic("secA2"[MM]),"_2:",sep=""), "r=0.042")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_MMAR_C2$model),
           col="red",
           lwd=2))

plot(data$Time, data$D2, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(atop(paste("Δ",italic("secA2"),"/p",italic("secA2"[MT]),"_2:",sep=""), "r=0.044")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_Rv_D2$model),
           col="red",
           lwd=2))

plot(data$Time, data$E2, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(paste("Δ",italic("esxBA"),"_2: r=0.041",sep="")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_ΔesxBA_E2$model),
           col="red",
           lwd=2))

plot(data$Time, data$A3, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(paste("WT_3: r=0.044",sep="")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_WT_A3$model),
           col="red",
           lwd=2))

plot(data$Time, data$B3, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(paste("Δ",italic("secA2"),"_3: r=0.031",sep="")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_ΔsecA2_B3$model),
           col="red",
           lwd=2))

plot(data$Time, data$C3, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(atop(paste("Δ",italic("secA2"),"/p",italic("secA2"[MM]),"_3:",sep=""), "r=0.044")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_MMAR_C3$model),
           col="red",
           lwd=2))

plot(data$Time, data$D3, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(atop(paste("Δ",italic("secA2"),"/p",italic("secA2"[MT]),"_3:",sep=""), "r=0.045")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_Rv_D3$model),
           col="red",
           lwd=2))

plot(data$Time, data$E3, pch=16, cex = 1, col = "black", xlab="Time (hours)", ylab="log Colony forming unit (CFU/mL)", 
     main=expression(paste("Δ",italic("esxBA"),"_3: r=0.043",sep="")), ylim=c(0,12),
     lines(data$Time,
           predict(gc_fit_ΔesxBA_E3$model),
           col="red",
           lwd=2))


#Add text as labels to the x and y axes of the graph
mtext("Time (hours)", side=1, line= 0, adj=0.5, outer=TRUE)
mtext("Optical Density (OD600)", line= 0, side=2, adj = 0.5, outer=TRUE)


dev.off()


par(mfrow=c(1,1))




#Make an output file and write the summarized data for all growth curves to it

GC_output <- "./GC_output.csv" 

#Replace ΔsecA2(8) with ΔsecA2
colnames(data)[colnames(data) == "Time"] <- "time"

gc_out <- SummarizeGrowthByPlate(data)


write.csv2(gc_out, file= GC_output, quote = FALSE, row.names = FALSE)



















  