
# Function GET_FITTEDSINGLEDRUGS()
get_fittedsingledrugs <- function(drug_response){	
		#Estimate the starting parameters ic50
drug_response = as.data.frame(apply(drug_response,2,as.numeric))
	estimate_param <- tryCatch({ 
      drm(effect ~ conc, data = drug_response, fct = LL.4(fixed = c(NA,NA,NA,NA),names = c("SLOPE","MIN","MAX","IC50")),na.action=na.omit,control = drmc(errorm = FALSE))
      }, warning=function(w){
        drm(effect ~ conc, data = drug_response, fct = L.4(fixed = c(NA,NA,NA,NA), names = c("SLOPE","MIN","MAX","IC50")),na.action=na.omit)
            
      },error=function(e){
			  drm(effect ~ conc, data = drug_response, fct = L.4(fixed = c(NA,NA,NA,NA),na.action=na.omit,names = c("SLOPE","MIN","MAX","IC50")))
			})
fitted(estimate_param)
# PR(estimate_param,drug_response$conc)
}
c.get_fittedsingledrugs <- cmpfun(get_fittedsingledrugs)

# Function SDBASELINECOR() 
# for baseline correction using single drugs
SDbaselineCor<-function(plate.mat,conc.range,pair.list){

conc.range = apply(conc.range,2,as.character)
pair.list$index = as.character(pair.list$index)
pm<-as.matrix(plate.mat)
pm=apply(pm,2,as.numeric)
pm[nrow(plate.mat),ncol(plate.mat)]<-NA

drugpair<-subset(pair.list,pair.list$index==1)
#print(drugpair)
d1c<-subset(conc.range,conc.range[,1]==drugpair$drug1)[,2:ncol(conc.range)]
d2c<-subset(conc.range,conc.range[,1]==drugpair$drug2)[,2:ncol(conc.range)] #conc test range for drug 2 (2nd col in vv/pairslist)
	
coleffect<-as.vector(as.matrix(pm)) #combined column wise (col1,col2,col3 etc)
roweffect<-as.vector(t(as.matrix(pm))) 
rownames(pm)<-as.character(d2c)
colnames(pm)<-as.character(d1c)

# print(pm)
d1<-as.vector(unlist(rep(d1c,times=ncol(pm)))) 
d2<-as.vector(unlist(rep(d2c,times=ncol(pm)))) 

concd1<-as.vector(sapply(d1c,rep,times=nrow(pm)))
concd2<-as.vector(sapply(d2c,rep,times=nrow(pm)))

rowwise<-data.frame(as.numeric(d1),as.numeric(concd2),roweffect) # data combined row wise, rows are concatenated
colwise<-data.frame(as.numeric(d2), as.numeric(concd1), coleffect) # data combined column wise, columns are concatenated

colnames(rowwise)[c(1,2)]<-c("conc_d1", "conc_d2")
colnames(colwise)[c(1,2)]<-c("conc_d2", "conc_d1")

single1<-rowwise[1:nrow(plate.mat),c(1,3)] # single drug denoted as the first row in the original matrix
colnames(single1)<-c("conc","effect")
# print(single1)

single2<-colwise[1:ncol(plate.mat),c(1,3)] # single drug denoted as the first column in the original matrix
colnames(single2)<-c("conc","effect")

#fitting single drugs independently
sd1<-c.get_fittedsingledrugs(single1[-1,])
# print(sd1)
# print(single1)
sd2<-c.get_fittedsingledrugs(single2[-1,])
# print(sd2)
# print(single2)

baseline<-(min(as.numeric(sd1))+min(as.numeric(sd2)))/2

# CORRECT MATRIX BY A WEIGHTED CORRECTION FACTOR
pm_cor<-pm-((100-pm)/100*baseline)

rownames(pm_cor)<-as.character(d2c)
colnames(pm_cor)<-as.character(d1c)

output = list(pm, pm_cor,drugpair)
return(output) # return the adjusted matrix, and drug names
}
c.SDbaselineCor<-cmpfun(SDbaselineCor)

#####################################TWO WAY FITTING ###########################################################
# ----------------------------
# new functions defined here
# ----------------------------
twowayfitting=function(cor_matrix,drug_pair) {
#print(cor_matrix)
#print(drug_pair)

# Fitting single drugs using logistic functions
# NA values treated
# Drug 1 fitting (fitting the first row)
#drug1=as.data.frame(mat.or.vec(7,0)) # first row in the matrix
drug1=as.data.frame(mat.or.vec(nrow(cor_matrix)-1,0)) # first row in the matrix
drug1$dose=as.numeric(colnames(cor_matrix)[-1])
#print(drug1)
drug1$logconc=log10(drug1$dose)
drug1$inhibition=as.numeric(cor_matrix[1,-1])
# NA values now treated
drug1_model=drm(inhibition ~ logconc, data = drug1, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10,na.action=na.omit)
drug1$fitted_inhibition = PR(drug1_model,log10(drug1$dose))

# Drug 2 fitting (fitting the first column)
#drug2=as.data.frame(mat.or.vec(7,0)) # first column in the matrix
drug2=as.data.frame(mat.or.vec(ncol(cor_matrix)-1,0)) # first column in the matrix
drug2$dose=as.numeric(rownames(cor_matrix)[-1])
#print(drug2)
drug2$logconc=log10(drug2$dose)
drug2$inhibition=as.numeric(cor_matrix[-1,1])
drug2_model=drm(inhibition ~ logconc, data = drug2, fct = L.4(fixed = c(NA, NA, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10,na.action=na.omit)
drug2$fitted_inhibition = PR(drug2_model,log10(drug2$dose))

# Update the first row and first column
cor_matrix_2 = mat.or.vec(nrow(cor_matrix),ncol(cor_matrix))
colnames(cor_matrix_2) = colnames(cor_matrix)
rownames(cor_matrix_2) = rownames(cor_matrix)
cor_matrix_2[1,c(2:ncol(cor_matrix))] = drug1$fitted_inhibition
cor_matrix_2[c(2:nrow(cor_matrix)),1] = drug2$fitted_inhibition
#print(cor_matrix_2)

# Update the column2-column8
cor_matrix_3 = cor_matrix_2
for (i in 2:ncol(cor_matrix)){
		tmp = as.data.frame(mat.or.vec(nrow(cor_matrix)-1,0))
		tmp$dose = drug2$dose
		tmp$logconc = drug2$logconc
		tmp$inhibition = cor_matrix[c(2:nrow(cor_matrix)),i]
		tmp_min = cor_matrix_2[1,i]
		tmp_model=drm(inhibition ~ logconc, data = tmp, fct = L.4(fixed = c(NA, tmp_min, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10,na.action=na.omit)
          tmp$fitted_inhibition = PR(tmp_model,tmp$logconc)
          if(tmp$fitted_inhibition[nrow(cor_matrix)-1]<0) tmp$fitted_inhibition[nrow(cor_matrix)-1]=tmp_min
          cor_matrix_3[c(2:nrow(cor_matrix)),i] = tmp$fitted_inhibition
          #print(tmp)
}
#print(cor_matrix_3)

# Update the row2-row8
cor_matrix_4 = cor_matrix_2
for (i in 2:nrow(cor_matrix)){
		tmp = as.data.frame(mat.or.vec(nrow(cor_matrix)-1,0))
		tmp$dose = drug1$dose
		tmp$logconc = drug1$logconc
		tmp$inhibition = cor_matrix[i,c(2:ncol(cor_matrix))]
		tmp_min = cor_matrix_2[i,1]
		tmp_model=drm(inhibition ~ logconc, data = tmp, fct = L.4(fixed = c(NA, tmp_min, NA,NA), names = c("SLOPE","MIN","MAX","IC50")),logDose=10,na.action=na.omit)
          tmp$fitted_inhibition = PR(tmp_model,tmp$logconc)
          if(tmp$fitted_inhibition[ncol(cor_matrix)-1]<0) tmp$fitted_inhibition[ncol(cor_matrix)-1]=tmp_min
          cor_matrix_4[i,c(2:ncol(cor_matrix))] = tmp$fitted_inhibition
          #print(tmp)
 }
#print(cor_matrix_4)
 
# take average of cor_matrix_3 and cor_matrix_4 as cor_matrix_final
cor_matrix_final = (cor_matrix_3+cor_matrix_4)/2
#print(cor_matrix_final)

# make cor_matrix_bliss based on cor_matrix_2
cor_matrix_bliss = cor_matrix_2
for (i in 2:nrow(cor_matrix_2)){
	for (j in 2:ncol((cor_matrix_2))){
		cor_matrix_bliss[i,j] = cor_matrix_2[i,1]+cor_matrix_2[1,j]-cor_matrix_2[i,1]*cor_matrix_2[1,j]/100
	}
}
# negative and positive controls are removed
cor_matrix_final[1,1]=0
cor_matrix_bliss[1,1]=0 
#cor_matrix_bliss[8,8]=0
cor_matrix_final[nrow(cor_matrix_final),ncol(cor_matrix_final)]=ifelse(cor_matrix_final[nrow(cor_matrix_final),ncol(cor_matrix_final)]>100,100,cor_matrix_final[nrow(cor_matrix_final),ncol(cor_matrix_final)]) # cannot be over 100 for the estimation

diff_cor_matrix=(cor_matrix_final-cor_matrix_bliss)

diff_cor_matrix_new<-diff_cor_matrix[-1,] ### remove first row
diff_cor_matrix_new<-diff_cor_matrix[,-1] ### remove first col
sum_diff_cor_matrix<-sum(diff_cor_matrix_new,na.rm=TRUE)/(ncol(diff_cor_matrix_new)*nrow(diff_cor_matrix_new)) 

output = sum_diff_cor_matrix
return(output)
} #END of the twowayfitting function

c.twowayfitting<-cmpfun(twowayfitting)


