##############################################################

### source code for matrixeQTL analysis during subsampling

##############################################################

# delete df.txt if it exists and prepare new analysis of subsampling
if (file.exists("df.txt")) {file.remove("df.txt")}  
### Run separately for every tissue tissue = 'Liver'; subset dataset
write.table(sub.mouse.liver.expression.eqtl, file = "expression/liver.expr.txt", 
    sep = "\t", row.names = FALSE, quote = FALSE)
# subset liver snp expression data
write.table(sub.BXD.geno.SNP.eqtl.for.liver, file = "genotypes/liver.snps.txt", 
    sep = "\t", row.names = FALSE, quote = FALSE)
file.remove("df.txt")
### step 1
tissue <- "liver"
### Running Matrix eQTL ###
library("MatrixEQTL")
### Load genotype info
snps <- SlicedData$new()
snps$LoadFile(paste0("genotypes/", tissue, ".snps.txt"), skipRows = 1, 
    skipColumns = 1, sliceSize = 500)
### Load gene expression info
expr <- SlicedData$new()
expr$LoadFile(paste0("expression/", tissue, ".expr.txt"), skipRows = 1, 
    skipColumns = 1, sliceSize = 500)
### Load covariates
cvrt <- SlicedData$new()
# cvrt$LoadFile(paste0('covariates/',tissue,'.covariates.txt'),
# skipRows = 1, skipColumns = 1, sliceSize = 500); Load gene locations
geneloc <- read.table(paste0("2016-09-08 ", tissue, ".gene.loc.txt"), sep = "\t", 
    header = TRUE, stringsAsFactors = FALSE)
### Load SNP locations
snpsloc <- read.table(paste0("2016-09-08 BXD.geno.loc.eqtl.for.", tissue, 
    ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
options(MatrixEQTL.dont.preserve.gene.object = TRUE)
### Run Matrix eQTL
me <- Matrix_eQTL_main(snps = snps, gene = expr, cvrt = cvrt, output_file_name = "", 
    pvOutputThreshold = 0, useModel = modelLINEAR, errorCovariance = numeric(), 
    verbose = TRUE, output_file_name.cis = paste0("eQTL_results_AL_", tissue, 
        "_cis.txt"), pvOutputThreshold.cis = 1, snpspos = snpsloc, genepos = geneloc, 
    cisDist = 1e+06, pvalue.hist = FALSE, noFDRsaveMemory = TRUE)
### Save the number of degrees of freedom for each tissue
cat(file = "df.txt", tissue, "\t", me$param$dfFull, "\n", append = TRUE)
tissue <- "lung"
### Running Matrix eQTL ###
library("MatrixEQTL")
### Load genotype info
snps <- SlicedData$new()
snps$LoadFile(paste0("genotypes/", tissue, ".snps.txt"), skipRows = 1, 
    skipColumns = 1, sliceSize = 500)
### Load gene expression info
expr <- SlicedData$new()
expr$LoadFile(paste0("expression/", tissue, ".expr.txt"), skipRows = 1, 
    skipColumns = 1, sliceSize = 500)
### Load covariates
cvrt <- SlicedData$new()
# cvrt$LoadFile(paste0('covariates/',tissue,'.covariates.txt'),
# skipRows = 1, skipColumns = 1, sliceSize = 500); Load gene locations
geneloc <- read.table(paste0("2016-09-08 ", tissue, ".gene.loc.txt"), sep = "\t", 
    header = TRUE, stringsAsFactors = FALSE)
### Load SNP locations
snpsloc <- read.table(paste0("2016-09-08 BXD.geno.loc.eqtl.for.", tissue, 
    ".txt"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
options(MatrixEQTL.dont.preserve.gene.object = TRUE)
### Run Matrix eQTL
me <- Matrix_eQTL_main(snps = snps, gene = expr, cvrt = cvrt, output_file_name = "", 
    pvOutputThreshold = 0, useModel = modelLINEAR, errorCovariance = numeric(), 
    verbose = TRUE, output_file_name.cis = paste0("eQTL_results_AL_", tissue, 
        "_cis.txt"), pvOutputThreshold.cis = 1, snpspos = snpsloc, genepos = geneloc, 
    cisDist = 1e+06, pvalue.hist = FALSE, noFDRsaveMemory = TRUE)
### Save the number of degrees of freedom for each tissue
cat(file = "df.txt", tissue, "\t", me$param$dfFull, "\n", append = TRUE)
### step2 Read df.txt for the list of tissues and degrees of freedom of
### linear models
df <- read.table("df.txt", stringsAsFactors = FALSE)
names(df) <- c("tissue", "df")
show(df)
### List vector for storing Matrix eQTL results
big.list <- vector("list", nrow(df))
### Store gene and SNP names from the first tissue for matching with
### other tissues
genes <- NULL
snps <- NULL
### colClasses for faster reading of Matrix eQTL output
cc.file <- NA
### Loop over tissues
for (t1 in 1:nrow(df)) {
    ### Get tissue name
    tissue <- df$tissue[t1]
    ### Load Matrix eQTL output for the given tissue
    start.time <- proc.time()[3]
    tbl <- read.table(paste0("eQTL_results_AL_", tissue, "_cis.txt"), header = T, 
        stringsAsFactors = FALSE, colClasses = cc.file)
    end.time <- proc.time()[3]
    cat(tissue, "loaded in", end.time - start.time, "sec.", nrow(tbl), 
        "gene-SNP pairs.", "\n")
    ### set colClasses for faster loading of other results
    if (any(is.na(cc.file))) {
        cc.file <- sapply(tbl, class)
    }
    ### Set gene and SNP names for matching
    if (is.null(snps)) 
        snps <- unique(tbl$SNP)
    if (is.null(genes)) 
        genes <- unique(tbl$gene)
    ### Match gene and SNP names from Matrix eQTL output to 'snps' and
    ### 'genes'
    gpos <- match(tbl$gene, genes, nomatch = 0L)
    spos <- match(tbl$SNP, snps, nomatch = 0L)
    ### Assign each gene-SNP pair a unique id for later matching with other tissues
    id <- gpos + 2 * spos * length(genes)
    ### Transform t-statistics into correlations
    r <- tbl$t.stat/sqrt(df$df[t1] + tbl$t.stat^2)
    ### Record id's and correlations
    big.list[[t1]] <- list(id = id, r = r)
    ### A bit of clean up to reduce memory requirements
    rm(tbl, gpos, spos, r, id, tissue, start.time, end.time)
    gc()
}
rm(t1, cc.file)
### Find the set of gene-SNP pairs present in results for all tissues
keep <- rep(TRUE, length(big.list[[1]]$id))
for (t1 in 2:nrow(df)) {
    mch <- match(big.list[[1]]$id, big.list[[t1]]$id, nomatch = 0L)
    keep[mch == 0] <- FALSE
    cat(df$tissue[t1], ", overlap size", sum(keep), "\n")
}
final.ids <- big.list[[1]]$id[keep]
rm(keep, mch, t1)
### Create and fill in the matrix of z-scores Z-scores are calculated
### from correlations
big.matrix <- matrix(NA_real_, nrow = length(final.ids), ncol = nrow(df))
fisher.transform <- function(r) {
    0.5 * log((1 + r)/(1 - r))
}
for (t1 in 1:nrow(df)) {
    mch <- match(final.ids, big.list[[t1]]$id)
    big.matrix[, t1] <- fisher.transform(big.list[[t1]]$r[mch]) * sqrt(df$df[t1] - 
        1)
    cat(t1, "\n")
}
stopifnot(!any(is.na(big.matrix)))
rm(t1, mch)
### Save the big matrix
save(list = "big.matrix", file = "z-score.matrix.Rdata", compress = FALSE)
### Save gene names and SNP names for rows of big matrix
writeLines(text = genes[final.ids%%(length(genes) * 2)], con = "z-score.matrix.genes.txt")
writeLines(text = snps[final.ids%/%(length(genes) * 2)], con = "z-score.matrix.snps.txt")
### step3 Set estimation parameters
maxIterations <- 100
### Load big matrix of z-scores
load(file = "z-score.matrix.Rdata")
dim(big.matrix)
### Initialize parameters
{
    param <- list()
    ### K - the number of tissues
    K <- ncol(big.matrix)
    ### Delta - null covariance matrix across tissues
    param$Delta <- matrix(0.05, K, K)
    diag(param$Delta) <- 1
    ### Sigma - signal covariance matrix across tissues
    param$Sigma <- matrix(3, K, K) + diag(K)
    ### P - the vector of probabilities
    param$P <- rep(1/2^K, 2^K)
    ### Psubs - the vector of active eQTLs for each element of P
    Psubs <- vector("list", 2^K)
    for (i in 1:2^K) {
        a <- 2^((K - 1):0)
        b <- 2 * a
        Psubs[[i]] <- as.double(((i - 1)%%b) >= a)
    }
    rm(a, b, i)
    param$Psubs <- Psubs
    rm(Psubs)
    ### loglik - the initial likelihood
    param$loglik <- -Inf
    rm(K)
}
### The function does a single iteration of the estimation procedure
DoIteration <- function(big.matrix, param) {
    ### extract current model parameters
    K <- ncol(big.matrix)
    m <- nrow(big.matrix)
    Delta <- param$Delta
    Sigma <- param$Sigma
    P <- param$P
    Psubs <- param$Psubs
    ### The function for matrix power
    mat.power <- function(mat, pow) {
        e <- eigen(mat)
        V <- e$vectors
        return(V %*% diag(e$values^pow) %*% t(V))
    }
    ### Start the timer
    tic <- proc.time()
    ### variables to accumulate loglik - likelihood newP - marginal
    ### probabilities newDelta - the new Delta matrix newSigmaPlusDelta - Delta+Sigma
    cum.loglik <- 0
    cum.newP <- 0
    cum.newDelta <- 0
    cum.newSigmaPlusDelta <- 0
    ### Do calculations in slices of 10000 gene-SNP pairs
    step1 <- 100000L
    for (j in 1:ceiling(m/step1)) {
        fr <- step1 * (j - 1) + 1
        to <- min(step1 * j, m)
        X <- big.matrix[fr:to, , drop = FALSE]
        ### likelihood for the slice
        prob <- matrix(0, nrow(X), length(P))
        for (i in 1:length(Psubs)) {
            sigma_star <- Delta + Sigma * tcrossprod(Psubs[[i]])
            sigma_hfiv <- mat.power(sigma_star, -0.5)
            sigma_dethfiv <- (det(sigma_star))^(-0.5)
            w <- (1/(2 * pi)^(K/2)) * (P[i] * sigma_dethfiv)
            prob[, i] <- exp(log(w) - colSums(tcrossprod(sigma_hfiv/sqrt(2), 
                X)^2))
        }
        cum.loglik <- cum.loglik + sum(log(rowSums(prob)))
        ### Normalize probabilities for each gene-SNP pair to add up to 1
        prob <- prob/rowSums(prob)
        ### new vector of P - tissue specificity probabilities
        cum.newP <- cum.newP + colSums(prob)
        cum.newDelta <- cum.newDelta + crossprod(X * sqrt(prob[, 1]))
        cum.newSigmaPlusDelta <- cum.newSigmaPlusDelta + crossprod(X * 
            sqrt(prob[, length(P)]))
    }
    {
        ### Calculate Delta from the cumulative sum
        Delta <- cum.newDelta/cum.newP[1]
        ### normalize to force the diagonal to 1
        Delta <- Delta * tcrossprod(sqrt(1/diag(Delta)))
        ### Same with Sigma
        Sigma <- cum.newSigmaPlusDelta/tail(cum.newP, 1) - Delta
        e <- eigen(Sigma)
        if (any(e$values < 0)) {
            Sigma <- e$vectors %*% diag(pmax(e$values, 0)) %*% t(e$vectors)
        }
    }
    P <- cum.newP/sum(cum.newP)
    toc <- proc.time()
    return(list(Delta = Delta, Sigma = Sigma, P = P, Psubs = Psubs, loglik = cum.loglik, 
        time = toc - tic))
}
### The 'paralist' list vector will store model estimates at each iteration
paralist <- vector("list", maxIterations + 1)
paralist[[1]] <- param
rm(param)
### Perform up to 'maxIterations' iteration
for (i in 2:length(paralist)) {
    paralist[[i]] <- DoIteration(big.matrix = big.matrix, param = paralist[[i - 
        1]])
    cat(i, "\t", paralist[[i]]$loglik - paralist[[i - 1]]$loglik, "\t", 
        paralist[[i]]$time[3], "\n")
    if (i > 10) 
        if (paralist[[i]]$loglik < paralist[[i - 1]]$loglik) 
            break
}
paralist <- paralist[!sapply(paralist, is.null)]
### Save the results
save(list = "paralist", file = "paralist.Rdata")
### step4 Parameters
local.FDR.threshold <- 1
output.file.name <- "MT-eQTLs.txt"
### Load big matrix of z-scores
load(file = "z-score.matrix.Rdata")
dim(big.matrix)
### Load gene names and SNP names matching the rows of big.matrix
gnames <- readLines("z-score.matrix.genes.txt")
snames <- readLines("z-score.matrix.snps.txt")
### Load tissue names
df <- read.table("df.txt", stringsAsFactors = FALSE)
names(df) <- c("tissue", "df")
show(df)
### Load parameter estimates and pick the last one
load("paralist.Rdata")
param <- tail(paralist, 1)[[1]]
### Number of tissues
K <- ncol(big.matrix)
m <- nrow(big.matrix)
### The function for matrix power
mat.power <- function(mat, pow) {
    e <- eigen(mat)
    V <- e$vectors
    return(V %*% diag(e$values^pow) %*% t(V))
}
### Matrix of possible tissue specificity profiles
Pmat <- simplify2array(param$Psubs)
### Call eQTLs and save in a file
fid <- file(description = output.file.name, open = "wt")
writeLines(con = fid, paste0("SNP\tgene\t", paste0("isEQTL.", df$tissue, 
    collapse = "\t"), "\t", paste0("marginalP.", df$tissue, collapse = "\t")))
### Do calculations in slices of 10000 gene-SNP pairs
step1 <- 10000L
cumdump <- 0
for (j in 1:ceiling(nrow(big.matrix)/step1)) {
    fr <- step1 * (j - 1) + 1
    to <- min(step1 * j, nrow(big.matrix))
    X <- big.matrix[fr:to, , drop = FALSE]
    ### likelihood for the slice
    prob <- matrix(0, nrow(X), length(param$P))
    for (i in 1:length(param$Psubs)) {
        sigma_star <- param$Delta + param$Sigma * tcrossprod(param$Psubs[[i]])
        sigma_hfiv <- mat.power(sigma_star, -0.5)
        sigma_dethfiv <- (det(sigma_star))^(-0.5)
        w <- (1/(2 * pi)^(K/2)) * (param$P[i] * sigma_dethfiv)
        prob[, i] <- exp(log(w) - colSums(tcrossprod(sigma_hfiv/sqrt(2), 
            X)^2))
    }
    prob <- prob/rowSums(prob)
    ### Select tests with eQTLs significant at local.FDR.threshold level
    keep <- (prob[, 1] <= local.FDR.threshold)
    if (any(keep)) {
        marginalProb <- tcrossprod(prob[keep, , drop = FALSE], 1 - Pmat)
        tissueSpecificity <- t(Pmat)[apply(X = prob[keep, , drop = FALSE], 
            MARGIN = 1, FUN = which.max), ]
        dump <- data.frame(snames[(fr:to)[keep]], gnames[(fr:to)[keep]], 
            tissueSpecificity, marginalProb, row.names = NULL, check.rows = FALSE, 
            check.names = FALSE, stringsAsFactors = FALSE)
        write.table(dump, file = fid, quote = FALSE, sep = "\t", row.names = FALSE, 
            col.names = FALSE)
    }
    cumdump <- cumdump + sum(keep)
    cat("Slice", j, "of", ceiling(nrow(big.matrix)/step1), " eQTLs recorded:", 
        cumdump, "\n")
}
close(fid)
### step5
MTeQTLs <- read.table(file = "MT-eQTLs.txt", header = T)
liver.ASE.ensembl <- read.table(file = "liver.ASE.ensembl.txt", header = T)
mouse430aensembl_id <- read.table(file = "2015-12-07 mouse430aensembl_id.txt", 
    header = T)
MTeQTLs <- merge(MTeQTLs, mouse430aensembl_id, by.x = "gene", by.y = "probe_id")
library(data.table)
MTeQTLs.min <- data.table(MTeQTLs, key = c("ensembl_id", "marginalP.liver"))
MTeQTLs.min <- MTeQTLs.min[J(unique(ensembl_id)), mult = "first"]
MTeQTLs_ASE <- MTeQTLs.min
MTeQTLs_ASE$ASE[MTeQTLs_ASE$ensembl_id %in% liver.ASE.ensembl$ensembl_gene_id] <- 1  # 1: ASE; 0: non-ASE 
MTeQTLs_ASE$ASE[!MTeQTLs_ASE$ensembl_id %in% liver.ASE.ensembl$ensembl_gene_id] <- 0  # 1: ASE; 0: non-ASE 
MTeQTLs_ASE <- subset(MTeQTLs_ASE, select = c("ensembl_id", "marginalP.liver", 
    "ASE", "gene", "SNP"))
MTeQTLs_ASE_3c <- subset(MTeQTLs_ASE, select = c("ensembl_id", "marginalP.liver", 
    "ASE"))
write.table(MTeQTLs_ASE_3c, file = "MTeQTLs_ASE_3c.txt", sep = "\t", row.names = FALSE, 
    quote = FALSE)
# delete df.txt and prepare new analysis of subsampling
file.remove("df.txt")
