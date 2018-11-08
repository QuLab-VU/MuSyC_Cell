rm(list=ls())
graphics.off()

# Load libraries required
library(drc)
library(caTools)
library(kriging)
library(lattice)
library(reshape2)
library(plotrix)
library(compiler)
library(truncnorm)

# Set working directory
setwd("/home/xnmeyer/Documents/Lab/Repos/synergy_theory/zip")
# Source the functions file 
source(paste(getwd(),"/functions.R",sep=""))

# Set the output folder 
cdir<-getwd()
myFolder<-c("Delta_score") ### if it is NULL, it automatically generates
if(!is.null(myFolder)){
	dir.create(file.path(cdir,myFolder))
	setwd(file.path(cdir,myFolder))
}

# Read response and metadata
response<-read.csv("../responses.csv")
metadata<-read.csv("../metadata.csv")

# get unique bolck ids
blockId<-unique(response$BlockId)

# Make matrix to store delta score
delta_score<-matrix(NA,length(blockId),4)
colnames(delta_score)<-c("blockId","Drug1","Drug2","Delta")

# Loop to compute delta for each combination
for(i in 1:length(blockId)){ #

blockIndex<-which(response$BlockId %in% blockId[i])
dataM<-response[blockIndex,c(2,3,4)]

#################### Preparing Input Files start #######################
#  Making plate matrix file
plate.file<-acast(dataM,Col~Row,value.var="Value")
plate.file<-apply(t(apply(t(plate.file),2,rev)),2,rev)
plate.mat <- 100-plate.file
plate.mat<-apply(plate.mat,2,as.numeric)


# Mapping concentration range file
meta.blockIndex<-which(metadata$BlockId %in% blockId[i])

drug1name<-as.character(as.matrix(metadata$RowName[meta.blockIndex])) ##colNames
drug2name<-"Ibrutinib" ###rowNames
colnames(plate.file)<-rep(drug1name,ncol(plate.file))
rownames(plate.file)<-rep(drug2name,nrow(plate.file))

drug1conc<-as.character(as.matrix(metadata$RowConcs[meta.blockIndex]))
drug1conc<-rev(as.numeric(unlist(strsplit(drug1conc,","))))*2
drug2conc<-as.character(as.matrix(metadata$ColConcs[meta.blockIndex]))
drug2conc<-rev(as.numeric(unlist(strsplit(drug2conc,","))))*2

conc.file<-rbind(drug1conc,drug2conc)
conc.file<-cbind(rbind(drug1name,drug2name),conc.file)
conc.range <-conc.file
conc.range[1,2:ncol(conc.range)]<-as.character(round(as.numeric(conc.range[1,2:ncol(conc.range)]),digits=1))
conc.range[2,2:ncol(conc.range)]<-as.character(round(as.numeric(conc.range[2,2:ncol(conc.range)]),digits=1))
# Add col names and row names to plate mat
colnames(plate.mat)<-colnames(plate.file)
rownames(plate.mat)<-rownames(plate.file)

# Making drug combination file
pairs.file<-matrix(NA,1,5)
pairs.file[1,]<-c(1,drug1name,drug2name,"cellline",1)
pairs.file <- as.data.frame(pairs.file)
pair.list <- pairs.file
colnames(pair.list)<-c("index","drug1","drug2","cell.line","plate")

#################### Preparing Input Files End #######################

####################### Baseline correction #####################

	output_baseline= SDbaselineCor(plate.mat,conc.range,pair.list)

######################## Two Way Fitting #######################
# Single plate analysis

# raw_matrix=output_baseline[[1]] # raw matrix
 cor_matrix=output_baseline[[2]] # matrix after baseline correction
 drug_pair=output_baseline[[3]] # drug names

 output = twowayfitting(cor_matrix,drug_pair)

 delta_score[i,1] = blockId[i]
 delta_score[i,2] = unique(colnames(plate.mat))
 delta_score[i,3] = unique(rownames(plate.mat))
 delta_score[i,4] = output

}
write.csv(delta_score,"Delta_score.csv",row.names=F)


