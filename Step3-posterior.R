liver.mouse.eQTL.bayesian<-read.table(file="liver.mouse.eQTL.bayesian.txt")
head(liver.mouse.eQTL.bayesian)


# Caculate variance for beta.tieda by following Brian Kulis' lecture notes
# Invert Tau and V
Tau_invert<-diag.inverse(Tau)
V_invert<-diag.inverse(V)
PS_invert<-Tau_invert + V_invert

# S in Brian Kulis' lecture note:PS
PS <- diag.inverse(PS_invert)
# retrieve posterior variance
ps<-diag(PS)
range(ps)

# reshape posterior variance to long format
ps.long <- melt(ps)
head(ps.long)
# Caculate sd: square root on variance
ps.long$betas.tieda.se<-(ps.long$value)^0.5
# combine sd to the data.frame
liver.mouse.eQTL.bayesian<-cbind(liver.mouse.eQTL.bayesian,ps.long$betas.tieda.se)

# head(liver.mouse.eQTL.bayesian)
# rename betas.tieda.se
liver.mouse.eQTL.bayesian<-rename(liver.mouse.eQTL.bayesian, c("ps.long$betas.tieda.se"="betas.tieda.se", "liver.beta_se"="betas.hat.se"))

#liver.mouse.eQTL.bayesian<-subset(liver.mouse.eQTL.bayesian, select = c("ensembl_id", "betas.hat", "betas.hat.se", "betas.tieda",
#"betas.tieda.se","liver_pvalue", "abs_lung.beta", "neg_log_liver_pvalue", "neg_log_lung_pvalue"))

# caculate probability of betas.tieda below 0 based on betas.tieda and standard deviation
liver.mouse.eQTL.bayesian$p.below.0 <- pnorm(0,liver.mouse.eQTL.bayesian$betas.tieda, liver.mouse.eQTL.bayesian$betas.tieda.se)

head(liver.mouse.eQTL.bayesian)
dim(liver.mouse.eQTL.bayesian)
summary(liver.mouse.eQTL.bayesian$betas.tieda.se)

range(liver.mouse.eQTL.bayesian$p.below.0)
write.table(liver.mouse.eQTL.bayesian,file="liver.mouse.eQTL.bayesian with beta.txt")

