liver.mouse.eQTL.bayesian.tau <- read.table("liver.mouse.eQTL.bayesian.tau.txt")
# chi-square (Fisher, 1932, Lancaster, 1961)
Fcomb <- function(ps) {
    k <- length(ps)
    temp <- -2 * sum(log(ps))
    pchisq(temp, 2 * k, lower.tail = F)
}
# normal (Liptak, 1958, Stouffer 1949)
Ncomb <- function(ps) {
    k <- length(ps)
    z <- qnorm((1 - ps))
    Ts <- sum(z)/sqrt(k)  # sum(1-Phi^-1(1-p))/sqrt(k)
    pnorm(Ts, lower.tail = F)  #Same as 1-Phi
}
# META
metapval <- apply(cbind(liver.mouse.eQTL.bayesian.tau$lung_pvalue, liver.mouse.eQTL.bayesian.tau$liver_pvalue), 
    1, Ncomb)
liver.mouse.eQTL.bayesian.tau$metapval <- metapval
# Multiple posterior prob by 2 CHANGE?
# liver.mouse.eQTL.bayesian.tau$p.below.0 =
# 2*liver.mouse.eQTL.bayesian.tau$p.below.0 MT Method
mtresults <- read.table(paste0("MTeQTLs_ASE_3c_", sebsetn, "s.txt"), header = TRUE)
minmtresults <- sapply(liver.mouse.eQTL.bayesian.tau$ensembl_id, function(x) min(mtresults[mtresults$ensembl_id == as.character(x), "marginalP.liver"]))
newmtresults <- data.frame(liver.mouse.eQTL.bayesian.tau$ensembl_id, minmtresults, liver.mouse.eQTL.bayesian.tau$eqtl)
colnames(newmtresults) <- c("ensembl_id", "marginalp", "eqtl")
newresults <- liver.mouse.eQTL.bayesian.tau[, c("ensembl_id", "lung_pvalue", "liver_pvalue", "metapval", "p.below.0", "eqtl")]
library(pROC)
# ROC plotting
pdf(paste0("subsampleproc1", sebsetn, ".pdf"), width = 4, height = 4)
rocobj1 <- plot.roc(newresults$eqtl, newresults$liver_pvalue, col = "black", legacy.axes = TRUE, yaxs = "i")
rocobj2 <- lines.roc(newresults$eqtl, newresults$p.below.0, col = "red")
rocobj3 <- lines.roc(newmtresults$eqtl, newmtresults$marginalp, col = "green")
rocobj4 <- lines.roc(newresults$eqtl, newresults$metapval, col = "blue")
rocobj5 <- lines.roc(newresults$eqtl, newresults$lung_pvalue, col = "purple")
legend("bottomright", legend = c("Conventional liver", "TA-eQTL", "MT", "Meta", "Conventional lung"), col = c("black", "red", "green", "blue", "purple"), lwd = 2, cex = 0.75, bty = "n")
dev.off()

# ROC plotting
pdf(paste0("subsampleproc14", sebsetn, ".pdf"), width = 4, height = 4)
rocobj1 <- plot.roc(newresults$eqtl, newresults$liver_pvalue, col = "black", legacy.axes = TRUE, yaxs = "i")
rocobj2 <- lines.roc(newresults$eqtl, newresults$p.below.0, col = "red")
rocobj3 <- lines.roc(newmtresults$eqtl, newmtresults$marginalp, col = "green")
rocobj4 <- lines.roc(newresults$eqtl, newresults$metapval, col = "blue")
rocobj5 <- lines.roc(newresults$eqtl, newresults$lung_pvalue, col = "purple")
legend("bottomright", legend = c("Conventional liver", "TA-eQTL", "MT", "Meta", "Conventional lung"), col = c("black", "red", "green", "blue", "purple"), lwd = 2, cex = 1.5, bty = "n")
dev.off()

orig_auc1 <- as.numeric(ci(newresults$eqtl, newresults$liver_pvalue))
bayesian_auc1 <- as.numeric(ci(newresults$eqtl, newresults$p.below.0))
lung_auc1 <- as.numeric(ci(newresults$eqtl, newresults$lung_pvalue))
meta_auc1 <- as.numeric(ci(newresults$eqtl, newresults$metapval))
mt_auc1 <- as.numeric(ci(newmtresults$eqtl, newmtresults$marginalp))
auc1 <- rbind(orig_auc1, bayesian_auc1, mt_auc1, meta_auc1, lung_auc1)
colnames(auc1) <- c("lowerCI", "mean", "upperCI")
auc1 <- data.frame(auc1)
auc2 <- round(auc1[, ], 2)
auc2$CI <- paste(auc2$lowerCI, auc2$upperCI, sep = ", ")
auc2$CI <- paste("(", auc2$CI, ")", sep = "")
auc2$lowerCI <- NULL
auc2$upperCI <- NULL
colnames(auc2) <- c("AUC", "CI")
rownames(auc2) <- c("Conventional liver", "TA-eQTL", "MT", "Meta", "Conventional lung")
auctable <- xtable(auc2)
print.xtable(auctable, type = "latex", file = paste0("auc", sebsetn, ".tex"), latex.environments = "center")
auc3 <- auc1
rownames(auc3) <- c("Conventional liver", "TA-eQTL", "MT", "Meta", "Conventional lung")
auc3$methods <- factor(row.names(auc3))
positions <- c("Conventional liver", "TA-eQTL", "MT", "Meta", "Conventional lung")
aucfivemethods <- ggplot(auc3, aes(x = methods, y = mean)) + geom_bar(stat = "identity", 
    fill = c("black", "red", "green", "blue", "purple")) + xlab("Methods") + 
    ylab("Area under the curve (AUC)") + geom_errorbar(aes(ymin = lowerCI, 
    ymax = upperCI), width = 0.1) + scale_x_discrete(limits = positions) + 
    coord_cartesian(ylim = c(0.5, 0.9)) + theme(axis.title.y = element_text(size = rel(1.8), 
    angle = 90)) + theme(axis.title.x = element_text(size = rel(1.8), angle = 0)) + 
    theme(axis.text.x = element_text(face = "bold", size = 12))
pdf("aucfivemethods.pdf")
print(aucfivemethods)
dev.off()


# significant testing to compare two ROC curves
orig.roc <- roc(newresults$eqtl, newresults$liver_pvalue)
bayesian.roc <- roc(newresults$eqtl, newresults$p.below.0)
mt.roc <- roc(newmtresults$eqtl, newmtresults$marginalp)
meta.roc <- roc(newresults$eqtl, newresults$metapval)
roc.test(orig.roc, bayesian.roc)
roc.test(orig.roc, mt.roc)
roc.test(orig.roc, meta.roc)
roc.test(mt.roc, bayesian.roc)
