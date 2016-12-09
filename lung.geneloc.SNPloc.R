# Prepare lung gene location and snp location for eQTL analysis
rm(list = ls())
setwd("/Volumes/Transcend/Thesis_project/eQTL data/1.2.lung.gene.snp")
BXD.geno <- read.table(file = "BXD-3.geno.txt", header = T)
# recode SNP, 'B to 0, D to 1'
BXD.geno1 <- as.data.frame(sapply(BXD.geno, gsub, pattern = "B", replacement = "0"))
BXD.geno2 <- as.data.frame(sapply(BXD.geno1, gsub, pattern = "H", replacement = "NA"))
BXD.geno3 <- as.data.frame(sapply(BXD.geno2, gsub, pattern = "D", replacement = "1"))
BXD.geno4 <- as.data.frame(sapply(BXD.geno3, gsub, pattern = "U", replacement = "NA"))
# SNPloc.txt was downloaded from BioMart-Ensembl website website
SNPloc <- read.table(file = "SNPloc1.txt", header = T)
SNPloc <- SNPloc[!duplicated(SNPloc$SNP), ]
SNPlibrary <- unique(SNPloc$SNP)
BXD.geno5 <- BXD.geno4[BXD.geno4$Locus %in% SNPlibrary, ]
BXD.geno.loc <- SNPloc[SNPloc$SNP %in% BXD.geno5$Locus, ]
BXD.geno.loc <- BXD.geno.loc[order(BXD.geno.loc$SNP), ]
# load lung expression data
mouse.lung.expression <- read.csv("GN160_DataAvgAnnot.rev0614.csv", 
    na.strings = c("", "NA"), header = TRUE, )
BXD.geno.SNP.library <- colnames(BXD.geno5)
mouse.lung.expression.eqtl <- mouse.lung.expression[, which(colnames(mouse.lung.expression) %in% BXD.geno.SNP.library)]
mouse.lung.expression.eqtl$Chr <- NULL
mouse.lung.expression.eqtl$Mb <- NULL
# reorder stain column names
mouse.lung.expression.eqtl <- mouse.lung.expression.eqtl[, order(colnames(mouse.lung.expression.eqtl))]
mouse.lung.expression.eqtl$ProbeSet <- mouse.lung.expression$ProbeSet
mouse.lung.expression.eqtl <- mouse.lung.expression.eqtl[, c(46, 1:45)]
# select SNP on the strains which has gene expression data available
BXD.geno6 <- BXD.geno5[, which(colnames(BXD.geno5) %in% colnames(mouse.lung.expression.eqtl))]
BXD.geno6 <- BXD.geno6[, order(colnames(BXD.geno6))]
BXD.geno6$Locus <- BXD.geno5$Locus
BXD.geno.SNP.eqtl <- BXD.geno6[, c(46, 1:45)]
BXD.geno.SNP.eqtl <- BXD.geno.SNP.eqtl[order(BXD.geno.SNP.eqtl$Locus), ]
# write BXD SNP genotypes for eqtl analysis
write.table(BXD.geno.SNP.eqtl, file = "2016-09-08 BXD.geno.SNP.eqtl.for.lung.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# write BXD SNP location for eqtl analysis
write.table(BXD.geno.loc, file = "2016-09-08 BXD.geno.loc.eqtl.for.lung.txt", sep = "\t", row.names = FALSE, quote = FALSE)
mouse430 <- read.table(file = "Affy mouse4302.txt", header = T)
mouse430 <- mouse430[!duplicated(mouse430$probeset), ]
lung.probeset.position.library <- mouse430$probeset
# subset mouse lung expression data with known gene location
mouse.lung.expression.eqtl.position <- mouse.lung.expression.eqtl[mouse.lung.expression.eqtl$ProbeSet %in% lung.probeset.position.library, ]
mouse.lung.expression.eqtl.position <- mouse.lung.expression.eqtl.position[order(mouse.lung.expression.eqtl.position$ProbeSet), ]
mouse.lung.expression.eqtl.position <- mouse.lung.expression.eqtl.position[order(mouse.lung.expression.eqtl.position$ProbeSet), ]
# write mouse lung gene expression data for eqtl analsis
write.table(mouse.lung.expression.eqtl.position, file = "2016-09-08 mouse.lung.expression.eqtl.txt", sep = "\t", row.names = FALSE, quote = FALSE)
lung.gene.loc <- mouse430[mouse430$probeset %in% mouse.lung.expression.eqtl.position$ProbeSet, ]
lung.gene.loc <- lung.gene.loc[order(lung.gene.loc$probeset), ]
write.table(lung.gene.loc, file = "2016-09-08 lung.gene.loc.txt", sep = "\t", row.names = FALSE, quote = FALSE)
