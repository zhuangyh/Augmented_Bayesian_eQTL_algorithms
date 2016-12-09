# Prepare liver gene location and snp location for eQTL analysis
rm(list = ls())
setwd("/Volumes/Transcend/Thesis_project/eQTL data/1.1.liver.gene.snp")
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
BXD.geno.SNP <- BXD.geno5[, c(2, 5:97)]
BXD.geno.SNP <- BXD.geno.SNP[order(BXD.geno.SNP$Locus), ]
BXD.geno.loc <- SNPloc[SNPloc$SNP %in% BXD.geno.SNP$Locus, ]
BXD.geno.loc <- BXD.geno.loc[order(BXD.geno.loc$SNP), ]
# Reformat mouse liver gene expression for Matrix eqtl analysis
mouse.liver.expression <- read.table("GN373_GSE16780_UCLA_Hybrid_MDP_Liver_Affy_HT_M430A_Sep11_RMA_Z-Score_Average.txt", 
    comment.char = "#", header = TRUE, sep = "\t", )
mouse.liver.expression <- as.data.frame(sapply(mouse.liver.expression, 
    gsub, pattern = "_at_A", replacement = "_at"))
# Creat the strain library with known SNP
BXD.geno.SNP.library <- colnames(BXD.geno.SNP)
mouse.liver.expression.eqtl <- mouse.liver.expression[, which(colnames(mouse.liver.expression) %in% BXD.geno.SNP.library)]
# reorder stain column names
mouse.liver.expression.eqtl <- mouse.liver.expression.eqtl[, order(colnames(mouse.liver.expression.eqtl))]
mouse.liver.expression.eqtl$ProbeSet <- mouse.liver.expression$ProbeSet
mouse.liver.expression.eqtl <- mouse.liver.expression.eqtl[, c(31, 1:30)]
# creat strain library with liver expression data
BXD.geno.strain.library <- colnames(mouse.liver.expression.eqtl)
# select SNP on the strains which has gene expression data available
BXD.geno.SNP1 <- BXD.geno.SNP[, which(colnames(BXD.geno.SNP) %in% BXD.geno.strain.library)]
BXD.geno.SNP1 <- BXD.geno.SNP1[, order(colnames(BXD.geno.SNP1))]
BXD.geno.SNP1$Locus <- BXD.geno.SNP$Locus
BXD.geno.SNP.eqtl <- BXD.geno.SNP1[, c(31, 1:30)]
BXD.geno.SNP.eqtl <- BXD.geno.SNP.eqtl[order(BXD.geno.SNP.eqtl$Locus), ]
# write BXD SNP genotypes for eqtl analysis
write.table(BXD.geno.SNP.eqtl, file = "2016-09-08 BXD.geno.SNP.eqtl.for.liver.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# check dimensions to make sure they match
dim(BXD.geno.loc)
# write BXD SNP location for eqtl analysis
write.table(BXD.geno.loc, file = "2016-09-08 BXD.geno.loc.eqtl.for.liver.txt", sep = "\t", row.names = FALSE, quote = FALSE)
# Affy_moe430a.txt was downloaded from BioMart-Ensembl website
mouse430a <- read.table(file = "Affy_moe430a1.txt", header = T)
mouse430a <- mouse430a[!duplicated(mouse430a$probeset), ]
liver.probeset.position.library <- mouse430a$probeset
# subset mouse liver expression data with known gene location
mouse.liver.expression.eqtl.position <- mouse.liver.expression.eqtl[mouse.liver.expression.eqtl$ProbeSet %in% 
    liver.probeset.position.library, ]
mouse.liver.expression.eqtl.position <- mouse.liver.expression.eqtl.position[order(mouse.liver.expression.eqtl.position$ProbeSet), ]
# write mouse liver gene expression data for eqtl analsis
write.table(mouse.liver.expression.eqtl.position, file = "2016-09-08 mouse.liver.expression.eqtl.txt",sep = "\t", row.names = FALSE, quote = FALSE)
liver.gene.loc <- mouse430a[mouse430a$probeset %in% mouse.liver.expression.eqtl.position$ProbeSet, ]
liver.gene.loc <- liver.gene.loc[order(liver.gene.loc$probeset), ]
write.table(liver.gene.loc, file = "2016-09-08 liver.gene.loc.txt", sep = "\t", row.names = FALSE, quote = FALSE)
