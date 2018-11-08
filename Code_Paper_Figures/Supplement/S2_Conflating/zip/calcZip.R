#! /usr/bin/Rscript
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
source(paste(getwd(),"/zip/functions.R",sep=""))

#################### Read in CSV files from python ##############
plate.mat = as.matrix(read.csv('zip/plate_mat.csv',header=T))
conc.range = read.csv('zip/conc_range.csv',header=T)
pair.list = read.csv('zip/pair_list.csv',header=T)
pair.list <- pair.list[,!(colnames(pair.list) %in% c("X"))]
####################### Baseline correction #####################
drug1name = as.character(pair.list$drug1)
drug2name = as.character(pair.list$drug2)
colnames(plate.mat)<-rep(drug1name,ncol(plate.mat))
rownames(plate.mat)<-rep(drug2name,nrow(plate.mat))
drug1conc<-conc.range$drug1name
drug2conc<-conc.range$drug2name
conc.file<-rbind(drug1conc,drug2conc)
conc.range<-cbind(rbind(drug1name,drug2name),conc.file)
output_baseline= SDbaselineCor(plate.mat,conc.range,pair.list)
######################## Two Way Fitting #######################
# Single plate analysis
# raw_matrix=output_baseline[[1]] # raw matrix
cor_matrix=output_baseline[[2]] # matrix after baseline correction
drug_pair=output_baseline[[3]] # drug names
output = twowayfitting(cor_matrix,drug_pair)
write.csv(output,"zip/delta_score.csv",row.names=F)


