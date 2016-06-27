####### Bayesian Method

# load mouse lung cis eqtl result
lung.mouse.eQTL<-read.table(file="mouselung.cis.1M.eqtls.txt",  header=T)
# load mouse liver cis eqtl result
liver.mouse.eQTL<-read.table(file="sub.mouseliver.cis.1M.eqtls.txt",  header=T)

mouse4302ensembl_id<-read.table(file="2015-12-04 mouse4302ensembl_id.txt",  header=T)
mouse430aensembl_id<-read.table(file="2015-12-07 mouse430aensembl_id.txt",  header=T)
# Add ensemble id annoatation to the data
lung.mouse.eQTL<-merge(lung.mouse.eQTL, mouse4302ensembl_id, by.x = "gene", by.y="probe_id")
liver.mouse.eQTL<-merge(liver.mouse.eQTL, mouse430aensembl_id, by.x = "gene", by.y="probe_id")
head(lung.mouse.eQTL)
head(liver.mouse.eQTL)

library(data.table)
library(plyr)
# Select lung Gene-SNP pair with minimum P value
lung.mouse.eQTL.min <- data.table(lung.mouse.eQTL, key=c('ensembl_id', "pvalue"))
lung.mouse.eQTL.min<-lung.mouse.eQTL.min[J(unique(ensembl_id)),mult="first"]
lung.mouse.eQTL.min<-as.data.frame(lung.mouse.eQTL.min)

# Select liver Gene-SNP pair with minimum P value
liver.mouse.eQTL.min <- data.table(liver.mouse.eQTL, key=c('ensembl_id', "pvalue"))
liver.mouse.eQTL.min<-liver.mouse.eQTL.min[J(unique(ensembl_id)),mult="first"]
liver.mouse.eQTL.min<-as.data.frame(liver.mouse.eQTL.min)



lung.mouse.eQTL.min<-rename(lung.mouse.eQTL.min, c("pvalue"="lung_pvalue", "beta"="lung.beta", "beta_se"="lung.beta_se"))
liver.mouse.eQTL.min<-rename(liver.mouse.eQTL.min, c("pvalue"="liver_pvalue", "beta"="liver.beta", "beta_se"="liver.beta_se"))

head(lung.mouse.eQTL.min)
head(liver.mouse.eQTL.min)
tail(liver.mouse.eQTL.min)
dim(lung.mouse.eQTL.min)
dim(liver.mouse.eQTL.min)
# lung, liver eqtl with ensemble_id
merged.mouse.eQTL.min<-merge(lung.mouse.eQTL.min, liver.mouse.eQTL.min, by.x = "ensembl_id", by.y="ensembl_id")
head(merged.mouse.eQTL.min)
dim(merged.mouse.eQTL.min)
merged.mouse.eQTL.min<-data.frame(merged.mouse.eQTL.min)
merged.mouse.eQTL.min<-merged.mouse.eQTL.min[, c(1, 5, 7, 8, 12, 14, 15 )]
head(merged.mouse.eQTL.min)
write.table(merged.mouse.eQTL.min,file="mouse.liver.expression.min.txt", sep="\t", row.names=FALSE, quote=FALSE)

library(fBasics)
table12 <- basicStats(merged.mouse.eQTL.min[, c(2:7)])[c("Mean", "Stdev", "Median", "Minimum", "Maximum"),]
table12 <- data.frame(table12)
table1 <- table12[, c(5, 6, 4)]
table2 <- table12[, c(2, 3, 1)]
library(plyr)
table1 <- rename(table1, c("liver_pvalue"="P value", "liver.beta"="beta", "liver.beta_se"="beta_se"))
table2 <- rename(table2, c("lung_pvalue"="P value", "lung.beta"="beta", "lung.beta_se"="beta_se"))

library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
table1 <- xtable(table1)
table2 <- xtable(table2)


print.xtable(table1, type="latex", file="table1.tex", latex.environments = "center")
print.xtable(table2, type="latex", file="table2.tex", latex.environments = "center")


Pvalue<-c(0.4, 0.2, 0.1, 0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000005)

chisq.results <- matrix(0, nrow=length(Pvalue), ncol=5)
colnames(chisq.results)<-c("Pvalue_threshold","Pvalue_chisq.test","actual_shared", "expected_shared", "Foldchange" )

# Populate the said matrix
for (i in 1:length(Pvalue)) {
  chisq.results[i,1] <- Pvalue[i]
  chisq.results[i,2] <- chisq.test(table(merged.mouse.eQTL.min$lung_pvalue<Pvalue[i], merged.mouse.eQTL.min$liver_pvalue<Pvalue[i]),correct=T)$p.value
  chisq.results[i,3] <- table(merged.mouse.eQTL.min$lung_pvalue<Pvalue[i], merged.mouse.eQTL.min$liver_pvalue<Pvalue[i])[2,2]
  chisq.results[i,4] <- chisq.test(table(merged.mouse.eQTL.min$lung_pvalue<Pvalue[i], merged.mouse.eQTL.min$liver_pvalue<Pvalue[i]),correct=T)$expected[2,2]
  chisq.results[i,5] <- chisq.results[i,3]/chisq.results[i,4]
}
print(chisq.results)

chisq.results.df<-as.data.frame(chisq.results)

library(ggplot2)
chisqfc <- ggplot(chisq.results.df, aes(x=-log(Pvalue_threshold), y=Foldchange)) +geom_point()+geom_smooth(method=lm)
pdf("chisqfc.pdf")
print(chisqfc)
dev.off()

chisqfctab <- xtable(chisq.results.df)
print.xtable(chisqfctab, type="latex", file="chisqfctab.tex", latex.environments = "center")


####### START HERE
merged.mouse.eQTL.min<-read.table(file="mouse.liver.expression.min.txt",  header=T)

###KK exploratory code
#plot(-log(merged.mouse.eQTL.min$lung_pvalue,10), -log(merged.mouse.eQTL.min$liver_pvalue,10))
#lungs = -log(merged.mouse.eQTL.min$lung_pvalue,10)
#livers = -log(merged.mouse.eQTL.min$liver_pvalue,10)
#mean(lungs[livers>10]>5)
#mean(livers[lungs>10]>5)

###KK added - didn't have abs beta variables

merged.mouse.eQTL.min$abs_liver.beta = abs(merged.mouse.eQTL.min$liver.beta)
merged.mouse.eQTL.min$abs_lung.beta = abs(merged.mouse.eQTL.min$lung.beta)
merged.mouse.eQTL.min$abs_liver.beta = abs(merged.mouse.eQTL.min$liver.beta)
merged.mouse.eQTL.min$abs_lung.beta = abs(merged.mouse.eQTL.min$lung.beta)
merged.mouse.eQTL.min$neg_log_lung_pvalue = -log10(merged.mouse.eQTL.min$lung_pvalue)
merged.mouse.eQTL.min$neg_log_liver_pvalue = -log10(merged.mouse.eQTL.min$liver_pvalue)

# Simple linear regression between abs_liver.beta and abs_lung.beta
# fit1<-summary(lm(abs_liver.beta ~ abs_lung.beta, data=merged.mouse.eQTL.min))
# fit1
# tau<-fit1$sigma**2
# check association between abs_liver.beta and abs.lung.beta

#Plots
#ggplot(merged.mouse.eQTL.min, aes(x=abs_lung.beta, y=abs_liver.beta)) +geom_point()+geom_smooth(method=lm)
#cor(merged.mouse.eQTL.min$abs_lung.beta, merged.mouse.eQTL.min$abs_liver.beta)
#ggplot(merged.mouse.eQTL.min, aes(x=lung.beta, y=liver.beta)) +geom_point()+geom_smooth(method=lm)
#cor(merged.mouse.eQTL.min$lung.beta, merged.mouse.eQTL.min$liver.beta)

merged.mouse.eQTL<-merged.mouse.eQTL.min
# retrieve ensembl_id
markers<-merged.mouse.eQTL[, 1]
# Yg=Ag + Bg*Xsnp+V
# retrieve betas.hat (liver.beta)
betas.hat<-merged.mouse.eQTL$abs_liver.beta
# retrieve liver.beta_se
se<-merged.mouse.eQTL$liver.beta_se

# create Z matrix with 2 columns: 1 for intercept,abs_lung.beta (merged.mouse.eQTL[,10])
Z<-as.matrix(merged.mouse.eQTL$abs_lung.beta)
Z<-as.matrix(merged.mouse.eQTL$neg_log_lung_pvalue) ##Use p-value as Z - didn't make a big difference
Z<-replace(Z,is.na(Z),0)
Z<-data.frame(1,Z) 
Z<-as.matrix(Z)
rowLength<-length(markers)

#CHANGE, include both beta and pvalue
#Z1<-as.matrix(merged.mouse.eQTL$abs_lung.beta)
#Z2<-as.matrix(merged.mouse.eQTL$neg_log_lung_pvalue)
#Z1<-replace(Z1, is.na(Z1),0)
#Z2<-replace(Z2, is.na(Z2),0)
#Z<-data.frame(1,Z1,Z2)
#Z<-as.matrix(Z)
#rowLength<-length(markers)

# Regression: abs_liver.beta = intercept + beta*abs_lung.beta + error
lmsummary<-summary(lm(abs_liver.beta~-1+Z, data=merged.mouse.eQTL))
lmsummary
model.prior = lm(abs_liver.beta~-1+Z, data=merged.mouse.eQTL)
# error ~ N(0, Tau)
tau<-lmsummary$sigma**2
tau
# output coeffieients (gamma matrix)
# gamma matrix
gamma<-as.matrix(lmsummary$coefficients[,1])
# transpose Z matrix
Z_transpose<-t(Z)
# create identity matrix
identity<-diag(nrow=rowLength)
# original betas.hat
betas.hat<-as.matrix(betas.hat)

#### WEIGHTS
useweights = 0 ##CHANGE TOGGLE
if(useweights ==1)
{	  
	  val = 1
	  weight = exp(-merged.mouse.eQTL.min$neg_log_lung_pvalue + val)
}

#create V matrix for liver_residual_variance
V <- matrix(0, rowLength, rowLength)
# V, liver residual variance
diag(V) <- merged.mouse.eQTL$liver.beta_se^2
# Creat Tau matrix
Tau<- diag(tau, rowLength, rowLength)
# follow Chen's paper and caculate s
s <-V + Tau
if(useweights ==1) {s <-V + diag(weight)*Tau}

# create inverse function for inversing diagnoal matrix
diag.inverse <- function(x){diag(1/diag(x), nrow(x), ncol(x))}
# create multiplication function for multiplicating two diagnoal matrix
diag.multi <- function(x,y){diag(diag(x)*diag(y), nrow(x), ncol(x))}
# inverse s
S <-diag.inverse(s)
# follow chen's paper to caculate omega
omega<-diag.multi(S, V)
# retrieve omega value from the matrix
omega.diag<-diag(omega )
# summary the omega value
summary(omega.diag)


#regression beta
regbeta <-Z %*% gamma
head(regbeta)
summary(regbeta)
# caculate betas.tieda with the formula in Chen's paper
constant = max(merged.mouse.eQTL.min$abs_liver.beta)/max(regbeta) ###CHANGE
betas.tieda<- constant * omega %*% Z %*% gamma + (identity-omega) %*% betas.hat
head(betas.tieda)
head(betas.hat)


markers1<-as.character(markers)
# combine ensemble_id, betas.hat and betas.tieda
outputVector<-c(markers1,betas.hat,betas.tieda)
write.table(matrix(outputVector,rowLength),file="hm_tau_hmresults.txt",col.names=FALSE,row.names=FALSE,quote=FALSE)
liver.mouse.eQTL.bayesian<-read.table(file="hm_tau_hmresults.txt")
colnames(liver.mouse.eQTL.bayesian)<-c( "ensembl_id", "betas.hat","betas.tieda")
head(liver.mouse.eQTL.bayesian)
# merge dataset with betas.hat and betas.tieda
liver.mouse.eQTL.bayesian<- merge(liver.mouse.eQTL.bayesian, merged.mouse.eQTL.min, by = "ensembl_id")
head(liver.mouse.eQTL.bayesian)

write.table(liver.mouse.eQTL.bayesian,file="liver.mouse.eQTL.bayesian.txt")

