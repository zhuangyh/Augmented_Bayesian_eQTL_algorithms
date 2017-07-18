##############################################################

### Step 3: Posterior estimation

##############################################################

liver.mouse.eQTL.bayesian <- read.table(file = "liver.mouse.eQTL.bayesian.txt")
# Caculate variance for beta.tieda by following Brian Kulis' lecture
# notes Invert Tau and V
Tau_invert <- diag.inverse(Tau)
V_invert <- diag.inverse(V)
PS_invert <- Tau_invert + V_invert
# S in Brian Kulis' lecture note:PS
PS <- diag.inverse(PS_invert)
# retrieve posterior variance
ps <- diag(PS)
# reshape posterior variance to long format
ps.long <- melt(ps)
# Caculate sd: square root on variance
ps.long$betas.tieda.se <- (ps.long$value)^0.5
# combine sd to the data.frame
liver.mouse.eQTL.bayesian <- cbind(liver.mouse.eQTL.bayesian, ps.long$betas.tieda.se)
# head(liver.mouse.eQTL.bayesian) rename betas.tieda.se
liver.mouse.eQTL.bayesian <- rename(liver.mouse.eQTL.bayesian, c(`ps.long$betas.tieda.se` = "betas.tieda.se", liver.beta_se = "betas.hat.se"))
# caculate probability of betas.tieda below 0 based on betas.tieda
# and standard deviation
liver.mouse.eQTL.bayesian$p.below.0 <- pnorm(0, liver.mouse.eQTL.bayesian$betas.tieda, liver.mouse.eQTL.bayesian$betas.tieda.se)
pdf("boxplotpb0.pdf")
boxplot(liver.mouse.eQTL.bayesian$p.below.0, ylab = "value", xlab = "Probability below 0")
dev.off()
write.table(liver.mouse.eQTL.bayesian, file = "liver.mouse.eQTL.bayesian with beta.txt")
Bayesianbetasd <- basicStats(liver.mouse.eQTL.bayesian[, c(3, 15, 16)])[c("Mean", "Stdev", "Median", "Minimum", "Maximum"), ]
Bayesianbetasd <- xtable(Bayesianbetasd)
print.xtable(Bayesianbetasd, type = "latex", file = "Bayesianbetasd.tex", 
    latex.environments = "center")
# Summary of posterior probability from weighted Bayesian method
eqtl.results1 <- matrix(0, nrow = length(Pthreshold), ncol = 3)
colnames(eqtl.results1) <- c("Pvalue_threshold", "Sig_in_liver", "Percent_in_liver")
for (i in 1:length(Pthreshold)) {
    eqtl.results1[i, 1] <- Pthreshold[i]
    eqtl.results1[i, 2] <- sum(liver.mouse.eQTL.bayesian$p.below.0 < 
        Pthreshold[i])
    eqtl.results1[i, 3] <- sum(liver.mouse.eQTL.bayesian$p.below.0 < 
        Pthreshold[i])/nrow(liver.mouse.eQTL.bayesian)
}
eqtl.results1 <- as.data.frame(eqtl.results1)
eqtl.results1$Pvalue_threshold <- as.character(eqtl.results1$Pvalue_threshold)
eqtl.results1$Sig_in_liver <- as.character(eqtl.results1$Sig_in_liver)
eqtl.results1$Percent_in_liver <- round(eqtl.results1$Percent_in_liver, 2)
weqtltab <- xtable(eqtl.results1)
print.xtable(weqtltab, type = "latex", include.rownames = FALSE, file = "weqtltab.tex", 
    latex.environments = "center")
