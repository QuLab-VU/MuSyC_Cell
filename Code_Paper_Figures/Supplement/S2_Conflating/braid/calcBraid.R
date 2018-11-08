#! /usr/bin/Rscript

library(braidrm)
df = read.csv('braid/effect.csv')
bd = braidrm(cbind(df$d1.conc,df$d2.conc),df$eff,fixed='ebraid')
write.csv(bd$coefficients[6],'braid/calc_kappa.csv')