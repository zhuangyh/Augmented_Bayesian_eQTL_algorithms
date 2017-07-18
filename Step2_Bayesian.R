##############################################################

### Step 2: Bayesian Method development
### Examine overlap of lung and liver cis-eQTL

##############################################################

lung.mouse.eQTL <- read.table(file = "mouselung.cis.1M.eqtls.txt", header = T)
# load mouse liver cis eqtl result
liver.mouse.eQTL <- read.table(file = "sub.mouseliver.cis.1M.eqtls.txt", header = T)
mouse4302ensembl_id <- read.table(file = "2015-12-04 mouse4302ensembl_id.txt", header = T)
mouse430aensembl_id <- read.table(file = "2015-12-07 mouse430aensembl_id.txt", header = T)
# Add ensemble id annoatation to the data
lung.mouse.eQTL <- merge(lung.mouse.eQTL, mouse4302ensembl_id, by.x = "gene", by.y = "probe_id")
liver.mouse.eQTL <- merge(liver.mouse.eQTL, mouse430aensembl_id, by.x = "gene", by.y = "probe_id")
library(data.table)
library(plyr)
# Select lung Gene-SNP pair with minimum P value
lung.mouse.eQTL.min <- data.table(lung.mouse.eQTL, key = c("ensembl_id", "pvalue"))
lung.mouse.eQTL.min <- lung.mouse.eQTL.min[J(unique(ensembl_id)), mult = "first"]
lung.mouse.eQTL.min <- as.data.frame(lung.mouse.eQTL.min)
# Select liver Gene-SNP pair with minimum P value
liver.mouse.eQTL.min <- data.table(liver.mouse.eQTL, key = c("ensembl_id", "pvalue"))
liver.mouse.eQTL.min <- liver.mouse.eQTL.min[J(unique(ensembl_id)), mult = "first"]
liver.mouse.eQTL.min <- as.data.frame(liver.mouse.eQTL.min)
lung.mouse.eQTL.min <- rename(lung.mouse.eQTL.min, c(pvalue = "lung_pvalue", beta = "lung.beta", beta_se = "lung.beta_se"))
liver.mouse.eQTL.min <- rename(liver.mouse.eQTL.min, c(pvalue = "liver_pvalue", beta = "liver.beta", beta_se = "liver.beta_se"))
# lung, liver eqtl with ensemble_id
merged.mouse.eQTL.min <- merge(lung.mouse.eQTL.min, liver.mouse.eQTL.min, by.x = "ensembl_id", by.y = "ensembl_id")
head(merged.mouse.eQTL.min)
dim(merged.mouse.eQTL.min)
merged.mouse.eQTL.min <- data.frame(merged.mouse.eQTL.min)
merged.mouse.eQTL.min <- merged.mouse.eQTL.min[, c(1, 5, 7, 8, 12, 14, 15)]
head(merged.mouse.eQTL.min)
write.table(merged.mouse.eQTL.min, file = "mouse.liver.expression.min.txt", 
    sep = "\t", row.names = FALSE, quote = FALSE)
Pthreshold <- c(0.05, 0.01, 0.001, 1e-04, 1e-05, 1e-06, 1e-07, 1e-08, 1e-09)
eqtl.results <- matrix(0, nrow = length(Pthreshold), ncol = 5)
colnames(eqtl.results) <- c("Pvalue_threshold", "Sig_in_lung", "Percent_in_lung", "Sig_in_liver", "Percent_in_liver")
# Populate the said matrix
for (i in 1:length(Pthreshold)) {
    eqtl.results[i, 1] <- Pthreshold[i]
    eqtl.results[i, 2] <- sum(merged.mouse.eQTL.min$lung_pvalue < Pthreshold[i])
    eqtl.results[i, 3] <- sum(merged.mouse.eQTL.min$lung_pvalue < Pthreshold[i])/nrow(merged.mouse.eQTL.min)
    eqtl.results[i, 4] <- sum(merged.mouse.eQTL.min$liver_pvalue < Pthreshold[i])
    eqtl.results[i, 5] <- sum(merged.mouse.eQTL.min$liver_pvalue < Pthreshold[i])/nrow(merged.mouse.eQTL.min)
}
eqtl.results <- as.data.frame(eqtl.results)
eqtl.results$Pvalue_threshold <- as.character(eqtl.results$Pvalue_threshold)
eqtl.results$Sig_in_lung <- as.character(eqtl.results$Sig_in_lung)
eqtl.results$Sig_in_liver <- as.character(eqtl.results$Sig_in_liver)
eqtl.results$Percent_in_lung <- round(eqtl.results$Percent_in_lung, 2)
eqtl.results$Percent_in_liver <- round(eqtl.results$Percent_in_liver, 2)
eqtltab <- xtable(eqtl.results)
print.xtable(eqtltab, type = "latex", include.rownames = FALSE, file = "eqtltab.tex", latex.environments = "center")
Pvalue <- c(seq(from = 0.95, to = 0.1, length.out = 18), 0.05, 0.01, 0.001, 1e-04, 1e-05, 1e-06, 1e-07, 5e-08)
chisq.results <- matrix(0, nrow = length(Pvalue), ncol = 5)
colnames(chisq.results) <- c("Pvalue_threshold", "Pvalue_chisq.test", "observed_shared", "expected_shared", "Folddifference")
# Populate the said matrix
for (i in 1:length(Pvalue)) {
    chisq.results[i, 1] <- Pvalue[i]
    a <- table(merged.mouse.eQTL.min$lung_pvalue<Pvalue[i], merged.mouse.eQTL.min$liver_pvalue<Pvalue[i])[,2]
    b <- chisq.test(table(merged.mouse.eQTL.min$lung_pvalue<Pvalue[i], merged.mouse.eQTL.min$liver_pvalue<Pvalue[i]),correct=T)$expected[, 2]
    c <- cbind(a,b)
    chisq.results[i,2] <- chisq.test(c,correct=T)$p.value
    chisq.results[i, 3] <- table(merged.mouse.eQTL.min$lung_pvalue < 
        Pvalue[i], merged.mouse.eQTL.min$liver_pvalue < Pvalue[i])[2, 2]
    chisq.results[i, 4] <- chisq.test(table(merged.mouse.eQTL.min$lung_pvalue < 
        Pvalue[i], merged.mouse.eQTL.min$liver_pvalue < Pvalue[i]), 
        correct = T)$expected[2, 2]
    chisq.results[i, 5] <- chisq.results[i, 3]/chisq.results[i, 4]
}
print(chisq.results)
chisq.results.df <- as.data.frame(chisq.results)

##############################################################

### Effect size and P value of cis-eQTLs in mouse lung and liver

##############################################################

library(ggplot2)
ae <- chisq.results.df[, c(1, 3, 4)]
ae1 <- data.frame(melt(ae, id.vars = "Pvalue_threshold"))
ae1$Pvalue_threshold <- as.numeric(ae1$Pvalue_threshold)
actvsexp <- ggplot(ae1, aes(x = -log10(Pvalue_threshold), y = value, color = variable)) + geom_line() + labs(y = "Number of overlapping cis-eQTL", 
    x = expression("-log"[10] ~ "(P value threshold)"))
# actvsexp1<- actvsexp + guides(fill=guide_legend(title=NULL))
actvsexp1 <- actvsexp + scale_colour_discrete(name = " ", breaks = c("observed_shared", 
    "expected_shared"), labels = c("observed overlap", "expected overlap")) + 
    scale_shape_discrete(name = " ", breaks = c("observed_shared", "expected_shared"), 
        labels = c("observed overlap", "expected overlap")) + geom_vline(xintercept = -log10(0.05), 
    color = "red", linetype = "dotted") + theme(legend.position = c(0.65, 0.8), text = element_text(size=15))
chisqfc <- ggplot(chisq.results.df, aes(x = -log10(Pvalue_threshold), 
    y = Folddifference)) + geom_point() + labs(y = "Ratio of Observed vs.Expected ", 
    x = expression("-log"[10] ~ "(P value threshold)")) + geom_hline(yintercept = 1, 
    color = "red", linetype = "dotted") + theme(text = element_text(size=15))
# Multiple plot function ggplot objects can be passed in ..., or to
# plotlist (as a list of ggplot objects) - cols: Number of columns
# in layout - layout: A matrix specifying the layout. If present,
# 'cols' is ignored.  If the layout is something like
# matrix(c(1,2,3,3), nrow=2, byrow=TRUE), then plot 1 will go in the
# upper left, 2 will go in the upper right, and 3 will go all the
# way across the bottom.
multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL) {
    library(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots <- length(plots)
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel ncol: Number of columns of plots nrow: Number of
        # rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)), ncol = cols, 
            nrow = ceiling(numPlots/cols))
    }
    if (numPlots == 1) {
        print(plots[[1]])
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this
            # subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row, 
                layout.pos.col = matchidx$col))
        }
    }
}
pdf("actvsexp.pdf", width = 7, height = 4.5)
multiplot(actvsexp1, chisqfc, cols = 2)
dev.off()
chisq.results.df$Pvalue_threshold <- as.character(chisq.results.df$Pvalue_threshold)
chisqfctab <- xtable(chisq.results.df, digits = c(0, 0, 4, 0, 0, 2))
chisqfctab
print.xtable(chisqfctab, type = "latex", file = "chisqfctab.tex", latex.environments = "center")

merged.mouse.eQTL.min <- read.table(file = "mouse.liver.expression.min.txt", 
    header = T)
merged.mouse.eQTL.min$abs_liver.beta <- abs(merged.mouse.eQTL.min$liver.beta)
merged.mouse.eQTL.min$abs_lung.beta <- abs(merged.mouse.eQTL.min$lung.beta)
merged.mouse.eQTL.min$abs_liver.beta <- abs(merged.mouse.eQTL.min$liver.beta)
merged.mouse.eQTL.min$abs_lung.beta <- abs(merged.mouse.eQTL.min$lung.beta)
merged.mouse.eQTL.min$neg_log_lung_pvalue <- -log10(merged.mouse.eQTL.min$lung_pvalue)
merged.mouse.eQTL.min$neg_log_liver_pvalue <- -log10(merged.mouse.eQTL.min$liver_pvalue)
cor.test(merged.mouse.eQTL.min$abs_lung.beta, merged.mouse.eQTL.min$neg_log_lung_pvalue)
# Make a basic volcano plot
vocano1 <- with(merged.mouse.eQTL.min, plot(lung.beta, -log10(lung_pvalue), 
    xlab = expression(beta[lgm]), ylab = expression("-log"[10] ~ "(lung P value)"), 
    xlim = c(-4, 4), ylim = c(0, 40)))
vocano2 <- with(merged.mouse.eQTL.min, plot(liver.beta, -log10(liver_pvalue), 
    xlab = expression(beta[vgm]), ylab = expression("-log"[10] ~ "(liver P value)")))
pdf("volcano.pdf", width = 8, height = 6)
par(mfrow = c(1, 2))
par(mar=c(5,5,2,2))
with(merged.mouse.eQTL.min, plot(lung.beta, -log10(lung_pvalue), xlab = expression(beta[lgm]), 
    ylab = expression("-log"[10] ~ "(lung P value)"), cex.lab= 2.2, xlim = c(-4, 4), ylim = c(0, 40)))
with(merged.mouse.eQTL.min, plot(liver.beta, -log10(liver_pvalue), xlab = expression(beta[vgm]), 
    ylab = expression("-log"[10] ~ "(liver P value)"), cex.lab= 2.2,xlim = c(-4, 4), ylim = c(0, 40)))
dev.off()
cor(merged.mouse.eQTL.min$abs_lung.beta, merged.mouse.eQTL.min$neg_log_lung_pvalue)
ggplot(merged.mouse.eQTL.min, aes(x = abs_lung.beta, y = abs_liver.beta)) + geom_point() + 
    xlab(expression(abs(hat(beta)) ~ "of lung")) + ylab(expression(abs(hat(beta)) ~ "of liver")) + 
    theme(text = element_text(size = 20)) + geom_abline(intercept = 0, slope = 1, colour = "red")                                                                                                                                             tilde(beta))) + theme(text = element_text(size = 20)) + geom_abline(intercept = 0, slope = 1, colour = "red")
ggplot(merged.mouse.eQTL.min, aes(x = lung.beta, y = liver.beta)) + geom_point() + 
  xlab(expression(abs(hat(beta)) ~ "of lung")) + ylab(expression(abs(hat(beta)) ~ "of liver")) + 
  theme(text = element_text(size = 20)) + geom_abline(intercept = 0, slope = 1, colour = "red")


##############################################################

### Bayesian Method development
### Unweighted Bayesian model

##############################################################

####### Start here for Bayesian analysis
merged.mouse.eQTL <- merged.mouse.eQTL.min
# retrieve ensembl_id
markers <- merged.mouse.eQTL[, 1]
# Yg=Ag + Bg*Xsnp+V retrieve betas.hat (liver.beta)
betas.hat <- merged.mouse.eQTL$abs_liver.beta
# retrieve liver.beta_se
se <- merged.mouse.eQTL$liver.beta_se
# create Z matrix with 2 columns: 1 for intercept,abs_lung.beta
# (merged.mouse.eQTL[,10])
Z <- as.matrix(merged.mouse.eQTL$abs_lung.beta)
Z <- as.matrix(merged.mouse.eQTL$neg_log_lung_pvalue)  ##Use p-value as Z - didn't make a big difference
Z <- replace(Z, is.na(Z), 0)
Z <- data.frame(1, Z)
Z <- as.matrix(Z)
rowLength <- length(markers)
# Regression: abs_liver.beta = intercept + beta*abs_lung.beta + error
lmsummary <- summary(lm(abs_liver.beta ~ -1 + Z, data = merged.mouse.eQTL))
lmsummary
model.prior <- lm(abs_liver.beta ~ -1 + Z, data = merged.mouse.eQTL)
# error ~ N(0, Tau)
tau <- lmsummary$sigma^2
# output coeffieients (gamma matrix) gamma matrix
gamma <- as.matrix(lmsummary$coefficients[, 1])
# transpose Z matrix
Z_transpose <- t(Z)
# create identity matrix
identity <- diag(nrow = rowLength)
# original betas.hat
betas.hat <- as.matrix(betas.hat)

##############################################################

### Bayesian Method development
### Weighted Bayesian model

##############################################################

#### WEIGHTS
useweights <- 0  ##CHANGE TOGGLE
if (useweights == 1) {
    val <- 1
    weight <- exp(-merged.mouse.eQTL.min$neg_log_lung_pvalue + val)
}
# create V matrix for liver_residual_variance
V <- matrix(0, rowLength, rowLength)
# V, liver residual variance
diag(V) <- merged.mouse.eQTL$liver.beta_se^2
# Creat Tau matrix
Tau <- diag(tau, rowLength, rowLength)
# follow Chen's paper and caculate s
s <- V + Tau
if (useweights == 1) {
    s <- V + diag(weight) * Tau
}
# create inverse function for inversing diagnoal matrix
diag.inverse <- function(x) {
    diag(1/diag(x), nrow(x), ncol(x))
}
# create multiplication function for multiplicating two diagnoal
# matrix
diag.multi <- function(x, y) {
    diag(diag(x) * diag(y), nrow(x), ncol(x))
}
# inverse s
S <- diag.inverse(s)
# follow chen's paper to caculate omega
omega <- diag.multi(S, V)
# retrieve omega value from the matrix
omega.diag <- diag(omega)
# summary the omega value
summary(omega.diag)
# regression beta
regbeta <- Z %*% gamma
summary(regbeta)
betas.tieda0 <- omega %*% Z %*% gamma + (identity - omega) %*% betas.hat
markers1 <- as.character(markers)
# combine ensemble_id, betas.hat and betas.tieda
outputVector <- c(markers1, betas.hat, betas.tieda0, regbeta)
write.table(matrix(outputVector, rowLength), file = "hm_tau_hmresults0.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
liver.mouse.eQTL.bayesian <- read.table(file = "hm_tau_hmresults0.txt")
colnames(liver.mouse.eQTL.bayesian) <- c("ensembl_id", "betas.hat", "betas.tieda", "regbeta")
head(liver.mouse.eQTL.bayesian)
# merge dataset with betas.hat and betas.tieda
liver.mouse.eQTL.bayesian <- merge(liver.mouse.eQTL.bayesian, merged.mouse.eQTL.min, by = "ensembl_id")
head(liver.mouse.eQTL.bayesian)
write.table(liver.mouse.eQTL.bayesian, file = "liver.mouse.eQTL.bayesian0.txt")
# plotting for Comparison of beta and posterior estimations in
# unweighted Bayesian model
unweighted <- ggplot(liver.mouse.eQTL.bayesian, aes(x = betas.hat, y = betas.tieda)) + geom_point() + xlab(expression(abs(hat(beta)))) + ylab(expression("unweighted " ~ 
    tilde(beta))) + theme(text = element_text(size = 40)) + geom_abline(intercept = 0, slope = 1, colour = "red")
pdf("unweighted.pdf")
print(unweighted)
dev.off()
# caculate betas.tieda with the formula in Chen's paper
constant <- max(merged.mouse.eQTL.min$abs_liver.beta)/max(regbeta)  ###CHANGE
betas.tieda <- constant * omega %*% Z %*% gamma + (identity - omega) %*% betas.hat
markers1 <- as.character(markers)
# combine ensemble_id, betas.hat and betas.tieda
outputVector <- c(markers1, betas.hat, betas.tieda, regbeta)
write.table(matrix(outputVector, rowLength), file = "hm_tau_hmresults.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
liver.mouse.eQTL.bayesian <- read.table(file = "hm_tau_hmresults.txt")
colnames(liver.mouse.eQTL.bayesian) <- c("ensembl_id", "betas.hat", "betas.tieda", "regbeta")
head(liver.mouse.eQTL.bayesian)
# merge dataset with betas.hat and betas.tieda
liver.mouse.eQTL.bayesian <- merge(liver.mouse.eQTL.bayesian, merged.mouse.eQTL.min, by = "ensembl_id")
write.table(liver.mouse.eQTL.bayesian, file = "liver.mouse.eQTL.bayesian.txt")
# plotting for Comparison of beta and posterior estimations in
# weighted Bayesian model
weighted <- ggplot(liver.mouse.eQTL.bayesian, aes(x = betas.hat, y = betas.tieda)) + geom_point() + xlab(expression(abs(hat(beta)))) + ylab(expression("weighted " ~ 
    tilde(beta))) + theme(text = element_text(size = 40)) + geom_abline(intercept = 0, slope = 1, colour = "red")
pdf("weighted.pdf")
print(weighted)
dev.off()
pdf("betacompa.pdf")
multiplot(unweighted, weighted, cols = 2)
dev.off()
