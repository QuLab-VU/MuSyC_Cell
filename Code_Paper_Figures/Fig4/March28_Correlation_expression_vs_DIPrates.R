#========================================================================================================
# Correlation between DIP rates and NOX5 expression. 
# 1. Find correlation in CCLE gene expression (public database). 
# 2. Validate the expression with qPCR in a panel of melanoma cell lines and isogenic sublines.
#========================================================================================================

source("/home/xnmeyer/Documents/Lab/Repos/drug_synergy_data_Cell/Code_Paper_Figures/Fig4/summarySE.R")
library(ggplot2)
dd = "/home/xnmeyer/Documents/Lab/Repos/drug_synergy_data_Cell/Code_Paper_Figures/Fig4"
setwd(dd)
par(ps = 8, cex = .5, cex.axis = .5)
#========================================================================================================
# Plot the correlation between expression and DIP rates (doublings/hr) at 8uM PLX4720
#========================================================================================================
dataE = read.csv("../../Data/NOX5_PGC1a_expression_qPCR_melanoma_panel.csv")

#========================================================================================================
# Plot for NOX5 expression and DIP rates
#========================================================================================================
data = dataE[(order(dataE$Cell)),]
data$logN = log2(data$rN)
df = summarySE(data, measurevar = "logN", groupvars = c("Cell"))
df$norm = unique(data$norm)
df$rates = unique(data$rates)
plot(logN~rates, data=df, cex=1,xlim= c(-0.035, 0.035),  ylim=c(-6, 6),
     ylab = "NOX5 Rel. Expression", xlab = "DIP Rates")
with(df[,], text(logN~rates, labels = (df$Cell), pos = 4, cex=0.5))
res2 = cor.test(df$rates, df$logN, method = c("pearson"))
res2
ggplot(data = df, aes(x = rates, y = logN)) + geom_point() + #main graph
  geom_errorbar(aes(ymin = logN-se, ymax = logN+se)) + 
  xlim(-0.025, 0.035) + ylim(-4, 5) + xlab("DIP Rates") +
  ylab("NOX5 expression") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=8)) +
  geom_text(aes(label=df$Cell),hjust=0, vjust=0,size=2.5) +
  ggsave(paste0("Correlation_NOX5_melanoma_qPCR.pdf"), width=2, height=1)
  
#========================================================================================================
# Plot for PGC1a expression and DIP rates
#========================================================================================================
data = dataE[(order(dataE$Cell)),]
data$logP = log2(data$rP)
df = summarySE(data, measurevar = "logP", groupvars = c("Cell"))
df$norm = unique(data$norm)
df$rates = unique(data$rates)
plot(logP~rates, data=df, cex=1,xlim= c(-0.035, 0.035),  ylim=c(-2.5, 2.5),
     ylab = "PGC1a Rel. Expression", xlab = "DIP Rates")
with(df[,], text(logP~rates, labels = (df$Cell), pos = 4, cex=0.5))
res2 = cor.test(df$rates, df$logP, method = c("pearson"))
res2
ggplot(data = df, aes(x = rates, y = logP)) + geom_point() + #main graph
  geom_errorbar(aes(ymin = logP-se, ymax = logP+se)) + 
  xlim(-0.025, 0.035) + ylim(-2, 1) + xlab("DIP Rates") +
  ylab("PGC1a expression") + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), text=element_text(size=8)) +
  geom_text(aes(label=df$Cell),hjust=0, vjust=0,size=2.5) +
  ggsave(paste0("Correlation_PGC1a_melanoma_qPCR.pdf"), width=2, height=1)


#========================================================================================================
# Plot for CCLE expression and DIP rates
#========================================================================================================

td = read.csv("../../Data/CCLE_expression_selected_cells.csv")

pdf("CCLE_correlation_NOX5.pdf", height = 4, width = 4)
plot(log2(NOX5)~rates, data=td, cex=0.8, xlim=c(-0.04, 0.04), ylim=c(0,10), 
     ylab = "NOX5 expression", xlab = "DIP rates", pch = 20)
# with(td[,], text(log2(as.integer(NOX5))~rates, labels = (td$cell), pos = 4, cex=0.5))
# summary(lm(log2(NOX5)~rates, data=td))
# abline(lm(log2(NOX5)~rates, data=td), lty=2)
# res1 = cor.test(td$rates, log2(td$NOX5), method = c("pearson"))
# text(0.02, 0.5, paste0("Corr = ", round(res1$estimate, 3)), cex = 0.8)
# res1
dev.off()

#========================================================================================================
#========================================================================================================