#=====================================================================================================
# Single-cell-derived SKMEL5 sublines initial response overlaid with clonal responses and population level
# concentration used: 8uM PLX4720
#=====================================================================================================

source("/home/xnmeyer/Documents/Lab/Repos/drug_synergy_data_Cell/Code_Paper_Figures/Fig4/summarySE.R")
library(ggplot2)
dd = "/home/xnmeyer/Documents/Lab/Repos/drug_synergy_data_Cell/Code_Paper_Figures/Fig4"
setwd(dd)
par(ps = 12, cex = 1, cex.axis = 1)

#=====================================================================================================
# data
#=====================================================================================================
# Clonal fractional proliferation cFP data at 8uM
data <- read.csv("../../Data/SKMEL5_cFP_assay_8uMPLX.csv", header=T, sep=",")
# Population level data t 8uM
pdat <- read.csv("../../Data/SKMEL5_population_8uMPLX.csv", header=T, sep=",")
pdat = summarySE(pdat, measurevar = "nl2", groupvars = c("Time"))
# Subclones initial screen at 8uM
cdat <- read.csv("../../Data/SKMEL5_subclones_8uMPLX.csv", header=T, sep=",")

#=====================================================================================================
# Plot
#=====================================================================================================

pdf(paste0("SKMEL5 + Sublines_Initial_Screen.pdf"), width=4, height=5)
cl1 <- grey.colors(length(unique(data$cID)))
par(ps = 12, cex = 1, cex.axis = 1)
plot(nl2~Time, data=data, type="n", ylim=c(-3.,2), xlim=c(0,105), ylab="", xlab="")
for(w in 1:length(unique(data$cID))){
  temp <- subset(data, data$cID==w)
  lines(temp$Time, temp$nl2, col=cl1[w])
}
lines(pdat$Time, pdat$nl2, col="black", lwd=3)
cl1 <- rainbow(length(unique(cdat$cID)))
par(ps = 12, cex = 1, cex.axis = 1)
#plot(nl2~Time, data=s, type="n", ylim=c(-1.2,1.2), xlim=c(0,80), ylab="", xlab="")
# for(w in (unique(cdat$cID))){
#   temp <- subset(cdat, cdat$cID==w)
#   lines(temp$Time, temp$nl2, col=cl1[w], lwd=1)
# }
lines(cdat$Time[cdat$cID==10], cdat$nl2[cdat$cID==10], col="deepskyblue", lwd=3)
lines(cdat$Time[cdat$cID==7], cdat$nl2[cdat$cID==7], col="green", lwd=3)
lines(cdat$Time[cdat$cID==1], cdat$nl2[cdat$cID==1], col="red", lwd=3)
dev.off()

#=====================================================================================================
#=====================================================================================================
