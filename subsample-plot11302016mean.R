rm(list = ls())
gc()
# set directory
setwd("/Volumes/Transcend/Thesis_project/eQTL data")
library(pROC)
library("MatrixEQTL")
library(fBasics)
library(plyr)
library(xtable)
library(data.table)
library(biomaRt)
library(ggplot2)
library(lme4)
library(lsmeans)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")
# subset dataset
combined_auc <- NULL
pdf("subsample.pdf", width = 8, height = 12)
par(mfrow = c(3, 2))
par(cex.lab= 2, cex.main=2)
# seed library
seedlib <- c(45:50)
aorig_auc <- abayesian_auc <- alung_auc <- ameta_auc <- amt_auc <- NULL
# set sub-sampling options: 10, 15, 20, 25, 30 strains
sublib <- c(25, 20, 15, 10, 30)
# loop to subsampling analyses
for (z in 1:length(sublib)) {
    sebsetn <- sublib[z]
    # full liver dataset has 30 strains
    combinded.results <- NULL
    orig_auc <- matrix(0, length(seedlib), 3)
    colnames(orig_auc) <- c("sebsetn", "samplingseed", "auc")
    bayesian_auc <- lung_auc <- meta_auc <- mt_auc <- orig_auc
    p <- 1
    for (k in 1:length(seedlib)) {
        set.seed(seedlib[k])
        # subset liver gene expression dataset
        mouse.liver.expression.eqtl <- read.table(file = "2016-09-08 mouse.liver.expression.eqtl.txt", 
            header = T)
        sub.mouse.liver.expression.eqtl <- mouse.liver.expression.eqtl[, 
            c(1, sample(2:dim(mouse.liver.expression.eqtl)[2], sebsetn, 
                replace = FALSE))]
        write.table(sub.mouse.liver.expression.eqtl, file = "sub.mouse.liver.expression.eqtl.txt", 
            sep = "\t", row.names = FALSE, quote = FALSE)
        # subset liver snp expression data
        BXD.geno.SNP.eqtl.for.liver <- read.table(file = "2016-09-08 BXD.geno.SNP.eqtl.for.liver.txt", header = T)
        head(BXD.geno.SNP.eqtl.for.liver)
        dim(BXD.geno.SNP.eqtl.for.liver)
        set.seed(seedlib[k])
        sub.BXD.geno.SNP.eqtl.for.liver <- BXD.geno.SNP.eqtl.for.liver[, c(1, sample(2:dim(BXD.geno.SNP.eqtl.for.liver)[2], sebsetn, replace = FALSE))]
        head(sub.BXD.geno.SNP.eqtl.for.liver)
        dim(sub.BXD.geno.SNP.eqtl.for.liver)
        write.table(sub.BXD.geno.SNP.eqtl.for.liver, file = "sub.BXD.geno.SNP.eqtl.for.liver.txt", sep = "\t", row.names = FALSE, quote = FALSE)
        ### MT eqtl analysis
        source("2016-09-12mtsubsetanalysis.R")
        ### liver eqtl analysis
        base.dir <- "/Volumes/Transcend/Thesis_project/eQTL data"
        # Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
        useModel <- modelLINEAR
        # Genotype file name
        SNP_file_name <- paste(base.dir, "/sub.BXD.geno.SNP.eqtl.for.liver.txt", sep = "")
        snps_location_file_name <- paste(base.dir, "/2016-09-08 BXD.geno.loc.eqtl.for.liver.txt", sep = "")
        # Gene expression file name
        expression_file_name <- paste(base.dir, "/sub.mouse.liver.expression.eqtl.txt", sep = "")
        gene_location_file_name <- paste(base.dir, "/2016-09-08 liver.gene.loc.txt", sep = "")
        # Covariates file name Set to character() for no covariates
        covariates_file_name <- character()
        # Output file name
        output_file_name_cis <- tempfile()
        output_file_name_tra <- tempfile()
        # Only associations significant at this level will be saved
        pvOutputThreshold_cis <- 1
        pvOutputThreshold_tra <- 5e-15
        # Error covariance matrix Set to numeric() for identity.
        errorCovariance <- numeric()
        # errorCovariance = read.table('Sample_Data/errorCovariance.txt');
        # Distance for local gene-SNP pairs
        cisDist <- 1e+06
        ## Load genotype data
        snps <- SlicedData$new()
        snps$fileDelimiter <- "\t"
        snps$fileOmitCharacters <- "NA"
        snps$fileSkipRows <- 1
        snps$fileSkipColumns <- 1
        snps$fileSliceSize <- 2000
        snps$LoadFile(SNP_file_name)
        ## Load gene expression data
        gene <- SlicedData$new()
        gene$fileDelimiter <- "\t"
        gene$fileOmitCharacters <- "NA"
        gene$fileSkipRows <- 1
        gene$fileSkipColumns <- 1
        gene$fileSliceSize <- 2000
        gene$LoadFile(expression_file_name)
        ## Load covariates
        cvrt <- SlicedData$new()
        cvrt$fileDelimiter <- "\t"
        cvrt$fileOmitCharacters <- "NA"
        cvrt$fileSkipRows <- 1
        cvrt$fileSkipColumns <- 1
        if (length(covariates_file_name) > 0) {
            cvrt$LoadFile(covariates_file_name)
        }
        ## Run the analysis
        snpspos <- read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
        genepos <- read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)
        head(genepos)
        me <- Matrix_eQTL_main(snps = snps, gene = gene, output_file_name = output_file_name_tra, 
            pvOutputThreshold = pvOutputThreshold_tra, useModel = useModel, 
            errorCovariance = numeric(), verbose = TRUE, output_file_name.cis = output_file_name_cis, 
            pvOutputThreshold.cis = pvOutputThreshold_cis, snpspos = snpspos, 
            genepos = genepos, cisDist = cisDist, pvalue.hist = TRUE, min.pv.by.genesnp = FALSE, noFDRsaveMemory = FALSE)
        unlink(output_file_name_cis)
        ## Results:
        cat("Analysis done in:", me$time.in.sec, " seconds", "\n")
        cat("Detected local eQTLs:", "\n")
        cis.eqtls <- me$cis$eqtls
        head(cis.eqtls)
        dim(cis.eqtls)
        cis.eqtls$beta_se <- cis.eqtls$beta/cis.eqtls$statistic
        write.table(cis.eqtls, file = "sub.mouseliver.cis.1M.eqtls.txt", sep = "\t", row.names = FALSE, quote = FALSE)
        ### eqtl analysis for lung Settings Linear model to use, modelANOVA,
        ### modelLINEAR, or modelLINEAR_CROSS
        useModel <- modelLINEAR
        # Genotype file name
        SNP_file_name <- paste(base.dir, "/2016-09-08 BXD.geno.SNP.eqtl.for.lung.txt", sep = "")
        snps_location_file_name <- paste(base.dir, "/2016-09-08 BXD.geno.loc.eqtl.for.lung.txt", sep = "")
        # Gene expression file name
        expression_file_name <- paste(base.dir, "/2016-09-08 mouse.lung.expression.eqtl.txt", sep = "")
        gene_location_file_name <- paste(base.dir, "/2016-09-08 lung.gene.loc.txt", sep = "")
        # Covariates file name Set to character() for no covariates
        covariates_file_name <- character()
        # Output file name
        output_file_name_cis <- tempfile()
        output_file_name_tra <- tempfile()
        # Only associations significant at this level will be saved
        pvOutputThreshold_cis <- 1
        pvOutputThreshold_tra <- 5e-15
        # Error covariance matrix Set to numeric() for identity.
        errorCovariance <- numeric()
        # errorCovariance = read.table('Sample_Data/errorCovariance.txt');
        # Distance for local gene-SNP pairs
        cisDist <- 1e+06
        ## Load genotype data
        snps <- SlicedData$new()
        snps$fileDelimiter <- "\t"
        snps$fileOmitCharacters <- "NA"
        snps$fileSkipRows <- 1
        snps$fileSkipColumns <- 1
        snps$fileSliceSize <- 2000
        snps$LoadFile(SNP_file_name)
        ## Load gene expression data
        gene <- SlicedData$new()
        gene$fileDelimiter <- "\t"
        gene$fileOmitCharacters <- "NA"
        gene$fileSkipRows <- 1
        gene$fileSkipColumns <- 1
        gene$fileSliceSize <- 2000
        gene$LoadFile(expression_file_name)
        ## Load covariates
        cvrt <- SlicedData$new()
        cvrt$fileDelimiter <- "\t"
        cvrt$fileOmitCharacters <- "NA"
        cvrt$fileSkipRows <- 1
        cvrt$fileSkipColumns <- 1
        if (length(covariates_file_name) > 0) {
            cvrt$LoadFile(covariates_file_name)
        }
        ## Run the analysis
        snpspos <- read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
        genepos <- read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)
        head(genepos)
        me <- Matrix_eQTL_main(snps = snps, gene = gene, output_file_name = output_file_name_tra, 
            pvOutputThreshold = pvOutputThreshold_tra, useModel = useModel, 
            errorCovariance = numeric(), verbose = TRUE, output_file_name.cis = output_file_name_cis, 
            pvOutputThreshold.cis = pvOutputThreshold_cis, snpspos = snpspos, 
            genepos = genepos, cisDist = cisDist, pvalue.hist = TRUE, min.pv.by.genesnp = FALSE, 
            noFDRsaveMemory = FALSE)
        unlink(output_file_name_cis)
        ## Results:
        cis.eqtls <- me$cis$eqtls
        head(cis.eqtls)
        dim(cis.eqtls)
        cis.eqtls$beta_se <- cis.eqtls$beta/cis.eqtls$statistic
        write.table(cis.eqtls, file = "mouselung.cis.1M.eqtls.txt", sep = "\t", row.names = FALSE, quote = FALSE)
        ####### Bayesian Method load mouse lung cis eqtl result
        lung.mouse.eQTL <- read.table(file = "mouselung.cis.1M.eqtls.txt", header = T)
        # load mouse liver cis eqtl result
        liver.mouse.eQTL <- read.table(file = "sub.mouseliver.cis.1M.eqtls.txt", header = T)
        mouse4302ensembl_id <- read.table(file = "2015-12-04 mouse4302ensembl_id.txt", header = T)
        mouse430aensembl_id <- read.table(file = "2015-12-07 mouse430aensembl_id.txt", header = T)
        # Add ensemble id annoatation to the data
        lung.mouse.eQTL <- merge(lung.mouse.eQTL, mouse4302ensembl_id, by.x = "gene", by.y = "probe_id")
        liver.mouse.eQTL <- merge(liver.mouse.eQTL, mouse430aensembl_id, by.x = "gene", by.y = "probe_id")
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
        merged.mouse.eQTL.min <- data.frame(merged.mouse.eQTL.min)
        merged.mouse.eQTL.min <- merged.mouse.eQTL.min[, c(1, 5, 7, 8, 12, 14, 15)]
        head(merged.mouse.eQTL.min)
        write.table(merged.mouse.eQTL.min, file = "mouse.liver.expression.min.txt", sep = "\t", row.names = FALSE, quote = FALSE)
        ####### START HERE
        merged.mouse.eQTL.min <- read.table(file = "mouse.liver.expression.min.txt", header = T)
        merged.mouse.eQTL.min$abs_liver.beta <- abs(merged.mouse.eQTL.min$liver.beta)
        merged.mouse.eQTL.min$abs_lung.beta <- abs(merged.mouse.eQTL.min$lung.beta)
        merged.mouse.eQTL.min$abs_liver.beta <- abs(merged.mouse.eQTL.min$liver.beta)
        merged.mouse.eQTL.min$abs_lung.beta <- abs(merged.mouse.eQTL.min$lung.beta)
        merged.mouse.eQTL.min$neg_log_lung_pvalue <- -log10(merged.mouse.eQTL.min$lung_pvalue)
        merged.mouse.eQTL.min$neg_log_liver_pvalue <- -log10(merged.mouse.eQTL.min$liver_pvalue)
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
        tau
        # output coeffieients (gamma matrix) gamma matrix
        gamma <- as.matrix(lmsummary$coefficients[, 1])
        # transpose Z matrix
        Z_transpose <- t(Z)
        # create identity matrix
        identity <- diag(nrow = rowLength)
        # original betas.hat
        betas.hat <- as.matrix(betas.hat)
        #### WEIGHTS
        useweights <- 0  ##CHANGE TOGGLE
        if (useweights == 1) {
            val <- 1
            weight <- exp(-merged.mouse.eQTL.min$neg_log_lung_pvalue + 
                val)
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
        # create multiplication function for multiplicating two diagnoal matrix
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
        betas.tieda0 <- omega %*% Z %*% gamma + (identity - omega) %*% 
            betas.hat
        markers1 <- as.character(markers)
        # combine ensemble_id, betas.hat and betas.tieda
        outputVector <- c(markers1, betas.hat, betas.tieda0, regbeta)
        write.table(matrix(outputVector, rowLength), file = "hm_tau_hmresults0.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
        liver.mouse.eQTL.bayesian <- read.table(file = "hm_tau_hmresults0.txt")
        colnames(liver.mouse.eQTL.bayesian) <- c("ensembl_id", "betas.hat", "betas.tieda", "regbeta")
        head(liver.mouse.eQTL.bayesian)
        # merge dataset with betas.hat and betas.tieda
        liver.mouse.eQTL.bayesian <- merge(liver.mouse.eQTL.bayesian, merged.mouse.eQTL.min, by = "ensembl_id")
        write.table(liver.mouse.eQTL.bayesian, file = "liver.mouse.eQTL.bayesian0.txt")
        # caculate betas.tieda with the formula in Chen's paper
        constant <- max(merged.mouse.eQTL.min$abs_liver.beta)/max(regbeta)  ###CHANGE
        betas.tieda <- constant * omega %*% Z %*% gamma + (identity - omega) %*% betas.hat
        markers1 <- as.character(markers)
        # combine ensemble_id, betas.hat and betas.tieda
        outputVector <- c(markers1, betas.hat, betas.tieda, regbeta)
        write.table(matrix(outputVector, rowLength), file = "hm_tau_hmresults.txt", col.names = FALSE, row.names = FALSE, quote = FALSE)
        liver.mouse.eQTL.bayesian <- read.table(file = "hm_tau_hmresults.txt")
        colnames(liver.mouse.eQTL.bayesian) <- c("ensembl_id", "betas.hat", "betas.tieda", "regbeta")
        # merge dataset with betas.hat and betas.tieda
        liver.mouse.eQTL.bayesian <- merge(liver.mouse.eQTL.bayesian, merged.mouse.eQTL.min, by = "ensembl_id")
        write.table(liver.mouse.eQTL.bayesian, file = "liver.mouse.eQTL.bayesian.txt")
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
        # rename betas.tieda.se
        liver.mouse.eQTL.bayesian <- rename(liver.mouse.eQTL.bayesian, c(`ps.long$betas.tieda.se` = "betas.tieda.se", liver.beta_se = "betas.hat.se"))
        liver.mouse.eQTL.bayesian$p.below.0 <- pnorm(0, liver.mouse.eQTL.bayesian$betas.tieda, liver.mouse.eQTL.bayesian$betas.tieda.se)
        write.table(liver.mouse.eQTL.bayesian, file = "liver.mouse.eQTL.bayesian with beta.txt")
        ### START HERE
        liver.mouse.eQTL.bayesian <- read.table(file = "liver.mouse.eQTL.bayesian with beta.txt")
        liver.mouse.eQTL.bayesian.tau <- liver.mouse.eQTL.bayesian
        ### ASE
        liver.ASE <- read.csv(file = "ASE.genetics.113.153882-6.csv")
        # 440 unique gene ID
        length(unique(liver.ASE$geneID))
        # verify ASE table
        liver.ASE1 <- liver.ASE[which(liver.ASE$replicate == "M.CH. DxB and BxD"), ]
        sub.liver.ASE <- liver.ASE1
        summary(sub.liver.ASE$pvalBH.DxB7)
        # sub.liver.ASE <- sub.liver.ASE[ sub.liver.ASE$geneID %in%
        # names(table(sub.liver.ASE$geneID))[table(sub.liver.ASE$geneID) >1] ,
        # ] check the remain gene number after subsetting
        liver.ASE.symbol <- unique(sub.liver.ASE$geneID)
        # Annoate gene symbol with ensemble.ID
        mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
        liver.ASE.ensembl <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), 
            filters = "mgi_symbol", values = liver.ASE.symbol, mart = mouse)
        liver.ASE.ensembl <- unique(liver.ASE.ensembl)
        # delete liver ASE ensemble ID which are not in the
        # liver.mouse.eQTL.bayesian data frame
        liver.ASE.ensembl <- liver.ASE.ensembl[liver.ASE.ensembl$ensembl_gene_id %in% liver.mouse.eQTL.bayesian.tau$ensembl_id, ]
        liver.mouse.eQTL.bayesian.tau$eqtl[liver.mouse.eQTL.bayesian.tau$ensembl_id %in% liver.ASE.ensembl$ensembl_gene_id] <- 1
        liver.mouse.eQTL.bayesian.tau$eqtl[!liver.mouse.eQTL.bayesian.tau$ensembl_id %in% liver.ASE.ensembl$ensembl_gene_id] <- 0
        write.table(liver.mouse.eQTL.bayesian.tau, "liver.mouse.eQTL.bayesian.tau.txt")
        liver.mouse.eQTL.bayesian.tau$neg_log_liver_pvalue <- -log10(liver.mouse.eQTL.bayesian.tau$liver_pvalue)
        by(liver.mouse.eQTL.bayesian.tau[, c(1, 7, 9, 14)], liver.mouse.eQTL.bayesian.tau[, "eqtl"], summary)
        liver.mouse.eQTL.bayesian.tau$eqtl[liver.mouse.eQTL.bayesian.tau$ensembl_id %in% liver.ASE.ensembl$ensembl_gene_id] <- 1
        liver.mouse.eQTL.bayesian.tau$eqtl[!liver.mouse.eQTL.bayesian.tau$ensembl_id %in% liver.ASE.ensembl$ensembl_gene_id] <- 0
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
        metapval <- apply(cbind(liver.mouse.eQTL.bayesian.tau$lung_pvalue, liver.mouse.eQTL.bayesian.tau$liver_pvalue), 1, Ncomb)
        liver.mouse.eQTL.bayesian.tau$metapval <- metapval
        # MT Method
        mtresults <- read.table(paste0("MTeQTLs_ASE_3c.txt"), header = TRUE)
        # Select Gene-SNP pair with minimum P value
        minmtresults <- sapply(liver.mouse.eQTL.bayesian.tau$ensembl_id, 
            function(x) min(mtresults[mtresults$ensembl_id == as.character(x), "marginalP.liver"]))
        newmtresults <- data.frame(liver.mouse.eQTL.bayesian.tau$ensembl_id, minmtresults, liver.mouse.eQTL.bayesian.tau$eqtl)
        colnames(newmtresults) <- c("ensembl_id", "marginalp", "eqtl")
        # Merge MT results with the other results
        newresults <- liver.mouse.eQTL.bayesian.tau[, c("ensembl_id", "lung_pvalue", "liver_pvalue", "metapval", "p.below.0", "eqtl")]
        newresults$marginalp <- newmtresults$marginalp
        # Merge results of 6 times randomizations
        combinded.results <- rbind(combinded.results, newresults)
        orig_auc[p, 1] <- sebsetn
        orig_auc[p, 2] <- seedlib[k]
        orig_auc[p, 3] <- auc(newresults$eqtl, newresults$liver_pvalue)
        bayesian_auc[p, 1] <- sebsetn
        bayesian_auc[p, 2] <- seedlib[k]
        bayesian_auc[p, 3] <- auc(newresults$eqtl, newresults$p.below.0)
        lung_auc[p, 1] <- sebsetn
        lung_auc[p, 2] <- seedlib[k]
        lung_auc[p, 3] <- auc(newresults$eqtl, newresults$lung_pvalue)
        meta_auc[p, 1] <- sebsetn
        meta_auc[p, 2] <- seedlib[k]
        meta_auc[p, 3] <- auc(newresults$eqtl, newresults$metapval)
        mt_auc[p, 1] <- sebsetn
        mt_auc[p, 2] <- seedlib[k]
        mt_auc[p, 3] <- auc(newresults$eqtl, newresults$marginalp)
        p <- p + 1
    }
    # Calcaculate means for roc curve plotting
    mean.results <- ddply(combinded.results, .(ensembl_id), summarize, 
        lung_pvalue = mean(lung_pvalue), liver_pvalue = mean(liver_pvalue), 
        metapval = mean(metapval), p.below.0 = mean(p.below.0), marginalp = mean(marginalp))
    mean.results$eqtl <- newresults$eqtl
    # Combine subsampling result
    aorig_auc <- rbind(aorig_auc, orig_auc)
    abayesian_auc <- rbind(abayesian_auc, bayesian_auc)
    alung_auc <- rbind(alung_auc, lung_auc)
    ameta_auc <- rbind(ameta_auc, meta_auc)
    amt_auc <- rbind(amt_auc, mt_auc)
    # ROC plotting
    rocobj1 <- plot.roc(mean.results$eqtl, mean.results$liver_pvalue, main = paste0(sebsetn, 
        " strains"), percent = TRUE, col = "black", legacy.axes = TRUE, yaxs = "i")
    rocobj2 <- lines.roc(mean.results$eqtl, mean.results$p.below.0, percent = TRUE, col = "red")
    rocobj3 <- lines.roc(mean.results$eqtl, mean.results$marginalp, percent = TRUE, col = "green")
    # legend(45,30, legend=c('Original liver', 'Bayesian', 'MT'),
    # col=c('black', 'red', 'green'), lwd=2, cex = 0.85, bty = 'n')
    legend("bottomright", legend = c("Conventional liver", "TA-eQTL", "MT"), 
        col = c("black", "red", "green"), lwd = 2, cex = 2, bty = "n")
}
dev.off()
aorig_auc <- data.frame(aorig_auc)
abayesian_auc <- data.frame(abayesian_auc)
alung_auc <- data.frame(alung_auc)
ameta_auc <- data.frame(ameta_auc)
amt_auc <- data.frame(amt_auc)
aorig_auc$methods <- "Conventional liver"
abayesian_auc$methods <- "TA-eQTL"
alung_auc$methods <- "Conventional lung"
ameta_auc$methods <- "Meta"
amt_auc$methods <- "MT"
comauc <- rbind(aorig_auc, abayesian_auc, amt_auc, ameta_auc, alung_auc)
write.table(comauc, "comauc0929.txt")
# comauc <- read.table("comauc0929.txt")
comauc$methods <- factor(comauc$methods, levels = unique(as.character(comauc$methods)))
## Gives count, mean, standard deviation, standard error of the mean,
## and confidence interval (default 95%).  data: a data frame.
## measurevar: the name of a column that contains the variable to be
## summariezed groupvars: a vector containing names of columns that
## contain grouping variables na.rm: a boolean that indicates whether to
## ignore NA's conf.interval: the percent range of the confidence
## interval (default is 95%)
summarySE <- function(data = NULL, measurevar, groupvars = NULL, na.rm = FALSE, 
    conf.interval = 0.95, .drop = TRUE) {
    library(plyr)
    # New version of length which can handle NA's: if na.rm==T, don't count
    # them
    length2 <- function(x, na.rm = FALSE) {
        if (na.rm) 
            sum(!is.na(x)) else length(x)
    }
    # This does the summary. For each group's data frame, return a vector
    # with N, mean, and sd
    datac <- ddply(data, groupvars, .drop = .drop, .fun = function(xx, 
        col) {
        c(N = length2(xx[[col]], na.rm = na.rm), mean = mean(xx[[col]], 
            na.rm = na.rm), sd = sd(xx[[col]], na.rm = na.rm), min = min(xx[[col]], 
            na.rm = na.rm), max = max(xx[[col]], na.rm = na.rm))
    }, measurevar)
    # Rename the 'mean' column
    datac <- rename(datac, c(mean = measurevar))
    datac$se <- datac$sd/sqrt(datac$N)  # Calculate standard error of the mean
    # Confidence interval multiplier for standard error Calculate
    # t-statistic for confidence interval: e.g., if conf.interval is .95,
    # use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + 0.5, datac$N - 1)
    datac$ci <- datac$se * ciMult
    return(datac)
}
sumcomauc <- summarySE(comauc, measurevar = "auc", groupvars = c("methods", "sebsetn"))
# Use sebsetn as a factor rather than numeric
sumcomauc2 <- sumcomauc
sumcomauc2$sebsetn <- factor(sumcomauc2$sebsetn)
aucsubrange <- ggplot(sumcomauc2, aes(x = sebsetn, y = auc, fill = methods)) + 
    geom_bar(position = position_dodge(), stat = "identity") + geom_errorbar(aes(ymin = min, 
    ymax = max), width = 0.2, position = position_dodge(0.9)) + xlab("Number of strains for analysis") + 
    ylab("AUC") + coord_cartesian(ylim = c(0.5, 1)) + scale_fill_manual(values = c("black", "red", "green", "blue", "purple"))
comauc$samplingseed <- factor(comauc$samplingseed)
sub10 <- subset(comauc, sebsetn == 10)
sub15 <- subset(comauc, sebsetn == 15)
sub20 <- subset(comauc, sebsetn == 20)
sub25 <- subset(comauc, sebsetn == 25)
# Mixed model with ramdom effect on samplingseed
test10 <- lmer(auc ~ methods + (1 | samplingseed), data = sub10)
test15 <- lmer(auc ~ methods + (1 | samplingseed), data = sub15)
test20 <- lmer(auc ~ methods + (1 | samplingseed), data = sub20)
test25 <- lmer(auc ~ methods + (1 | samplingseed), data = sub25)
# Pair comparisons between methods
lsmeans(test10, pairwise ~ methods)
lsmeans(test15, pairwise ~ methods)
lsmeans(test20, pairwise ~ methods)
lsmeans(test25, pairwise ~ methods)
# Normalized AUC
comauc.liver <- subset(comauc, methods == "Conventional liver")
comauc.liver$liverauc <- comauc.liver$auc
comauc.liver$auc <- NULL
comauc.liver$methods <- NULL
ncomauc <- merge(comauc, comauc.liver, by = c("sebsetn", "samplingseed"))
ncomauc$nauc <- ncomauc$auc/ncomauc$liverauc
sumncomauc <- summarySE(ncomauc, measurevar = "nauc", groupvars = c("methods", 
    "sebsetn"))
# Use sebsetn as a factor rather than numeric
sumncomauc2 <- sumncomauc
sumncomauc2$sebsetn <- factor(sumncomauc2$sebsetn)
naucsubrange <- ggplot(sumncomauc2, aes(x = sebsetn, y = nauc, fill = methods)) + 
    geom_bar(position = position_dodge(), stat = "identity") + geom_errorbar(aes(ymin = min, 
    ymax = max), width = 0.2, position = position_dodge(0.9)) + xlab("Number of strains for analysis") + 
    ylab(expression("Ratio of AUC vs AUC"["Conventional liver"])) + coord_cartesian(ylim = c(0.8, 
    1.2)) + scale_fill_manual(values = c("black", "red", "green", "blue", "purple"))
pdf("naucsubrange.pdf", width = 8, height = 4)
print(naucsubrange)
dev.off()
sumncomauc2$min <- NULL
sumncomauc2$max <- NULL
sumcomauc2$min <- NULL
sumcomauc2$max <- NULL
sumcomauc_wide1 <- dcast(sumcomauc2, methods ~ sebsetn, value.var = "auc")
sumcomauc_wide2 <- dcast(sumncomauc2, methods ~ sebsetn, value.var = "nauc")
sumcomauc_wide <- cbind(sumcomauc_wide1, sumcomauc_wide2[, 2:5])
sumcomauc_wide <- sumcomauc_wide[, c(1, 2, 6, 3, 7, 4, 8, 5, 9)]
colnames(sumcomauc_wide) <- c("methods", "10_mean", "10_ratio", "15_mean", "15_ratio", "20_mean", "20_ratio", "25_mean", "25_ratio")
write.table(sumcomauc_wide, "sumcomauc_wide0929.txt")
print.xtable(xtable(sumcomauc_wide), type = "latex", file = "combined_auc.tex", 
    latex.environments = "center", include.rownames = FALSE)
