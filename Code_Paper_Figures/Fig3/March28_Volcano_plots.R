# =========================================================================================================
# Generate Volcano plot for differentially expressed genes between SKMEL5 subclones 
# SC01 vs SC10; SC01 vs SC07 and SC07 vs SC10
# =========================================================================================================
# Set the working directory
# =========================================================================================================

WORKDIR <- '/Users/paudelbb/RNA_Seq_Samples/Synergy_paper/data'
setwd(WORKDIR)

# =========================================================================================================
# Select the file for analysis and visualization
# =========================================================================================================
# Read DEGS files for analysis
filename = basename(file.choose())
res <- read.csv(filename, header=T)

# =========================================================================================================
# extract names for titles of the volcano plots
# =========================================================================================================
names = substr(filename, 6, 17)
head(res)
res = na.omit(res)
#res = subset(res, res$padj <= 0.01)
res$pvalue[res$pvalue==0] = 10^-100

# =========================================================================================================
# Make a basic volcano plot
# =========================================================================================================
value = 2
res1 = subset(res, res$padj <= 0.0001)
up = subset(res1, res1$log2FoldChange >=   value)
dn = subset(res1, res1$log2FoldChange <= - value)
set.seed(100)
upgenes = sample(up$symbol, 25, replace = F)
dngenes = sample(dn$symbol, 25, replace = F)
# =========================================================================================================
# =========================================================================================================

pdf(paste0(WORKDIR, "/",names, "_volcano_plot",".pdf"), width=5.5, height=5.5)
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.0, cex.axis=1.0, cex.lab=1.0)

with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main=paste0(names), xlim=c(-10,10), ylim=c(0, 350), 
               ylab="", col="black", cex=0.2))
with(subset(res, padj<= 0.001 & abs(log2FoldChange)>=value), points(log2FoldChange, -log10(pvalue), 
                                                                    pch=20, col="grey", cex=0.5))

# upregulated genes
sel <- res[(res$symbol %in% upgenes),] # or whatever you want to use
points(sel$log2FoldChange, -log10(sel$padj), col="brown", cex=0.3)
text(x=(sel$log2FoldChange) , y=(-log10(sel$pvalue)),
     label=sel$symbol, cex=0.5, col = "brown")
#downregulated genes
sel <- res[(res$symbol %in% dngenes),] # or whatever you want to use
points(sel$log2FoldChange, -log10(sel$padj),col="blue", cex=0.3)
text(x=(sel$log2FoldChange-1) , y=(-log10(sel$pvalue)),
     label=sel$symbol, cex=0.5, col = "blue")

#selected genes
sign.genes=c("NOX5", "PPARGC1A")
symbols(x=res$log2FoldChange[res$symbol == sign.genes[1]] , y=-log10(res$pvalue[res$symbol == sign.genes[1]]), 
        circles = c(0.01), inches = 0.04, bg="red", add=T)
text(x=res$log2FoldChange[res$symbol == sign.genes[1]] , y=(-log10(res$pvalue[res$symbol == sign.genes[1]])+15),
     label=sign.genes[1], cex=0.8, col = "red")

symbols(x=res$log2FoldChange[res$symbol == sign.genes[2]] , y=-log10(res$pvalue[res$symbol == sign.genes[2]]), 
        circles = c(0.01), inches = 0.04, bg="blue", add=T)
text(x=res$log2FoldChange[res$symbol == sign.genes[2]] , y=(-log10(res$pvalue[res$symbol == sign.genes[2]])-15), 
     label=sign.genes[2], cex=0.8, col = "blue")
abline(v=c(-value, value), lty=2)

# number of genes up or downregulated
text(-5 ,330,  label = paste0("n = ", dim(dn)[1]), cex=0.8, col = "blue")
text( 5 , 330, label = paste0("n = ", dim(up)[1]), cex=0.8, col = "brown")

dev.off()
# =========================================================================================================
# =========================================================================================================
