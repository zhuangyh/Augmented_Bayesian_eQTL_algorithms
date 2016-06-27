rm(list = ls())
gc()

#set directory
setwd("/Volumes/Transcend/Thesis_project/Subsetted_liver")
# subset dataset
sebsetn <- 30# full liver dataset has 30 strains

# subset liver gene expression dataset
mouse.liver.expression.eqtl <-read.table(file="2016-05-16 mouse.liver.expression.eqtl.txt",  header=T)
set.seed(50)
sub.mouse.liver.expression.eqtl <- mouse.liver.expression.eqtl[, c(1, sample(2:dim(mouse.liver.expression.eqtl)[2],sebsetn, replace=FALSE))]
write.table(sub.mouse.liver.expression.eqtl,file="sub.mouse.liver.expression.eqtl.txt", sep="\t", row.names=FALSE, quote=FALSE)

#subset liver snp expression data
BXD.geno.SNP.eqtl.for.liver <-read.table(file="2016-05-16 BXD.geno.SNP.eqtl.for.liver.txt",  header=T)
head(BXD.geno.SNP.eqtl.for.liver)
dim(BXD.geno.SNP.eqtl.for.liver)
set.seed(50)
sub.BXD.geno.SNP.eqtl.for.liver <- BXD.geno.SNP.eqtl.for.liver[, c(1, sample(2:dim(BXD.geno.SNP.eqtl.for.liver)[2],sebsetn, replace=FALSE))]
head(sub.BXD.geno.SNP.eqtl.for.liver)
dim(sub.BXD.geno.SNP.eqtl.for.liver)
write.table(sub.BXD.geno.SNP.eqtl.for.liver,file="sub.BXD.geno.SNP.eqtl.for.liver.txt", sep="\t", row.names=FALSE, quote=FALSE)

##########################

# liver eqtl analysis
library("MatrixEQTL")
base.dir = "/Volumes/Transcend/Thesis_project/Subsetted_liver"
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
# Genotype file name
SNP_file_name = paste(base.dir, "/sub.BXD.geno.SNP.eqtl.for.liver.txt", sep="");
snps_location_file_name = paste(base.dir, "/2016-05-16 BXD.geno.loc.eqtl.for.liver.txt", sep="");
# Gene expression file name
expression_file_name = paste(base.dir, "/sub.mouse.liver.expression.eqtl.txt", sep="");
gene_location_file_name = paste(base.dir, "/2016-05-16 liver.gene.loc.txt", sep="");
# Covariates file name
# Set to character() for no covariates
covariates_file_name = character() ;
# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0.000000000000005;
# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");
# Distance for local gene-SNP pairs
cisDist = 1e6; ##### 1 MB

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;
snps$fileSkipColumns = 1;
snps$fileSliceSize = 2000;
snps$LoadFile(SNP_file_name);

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;
gene$fileSkipColumns = 1;
gene$fileSliceSize = 2000;
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
cvrt$LoadFile(covariates_file_name);
}

## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
head(genepos)
me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_cis);
## Results:
cat('Analysis done in:', me$time.in.sec, ' seconds', '\n')
cat('Detected local eQTLs:','\n')
cis.eqtls<-me$cis$eqtls
head(cis.eqtls)
dim(cis.eqtls)
cis.eqtls$beta_se <-cis.eqtls$beta/cis.eqtls$statistic
write.table(cis.eqtls,file="sub.mouseliver.cis.1M.eqtls.txt", sep="\t", row.names=FALSE, quote=FALSE)

#########################
# eqtl analysis for lung
## Settings
# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
# Genotype file name
SNP_file_name = paste(base.dir, "/2016-05-16 BXD.geno.SNP.eqtl.for.lung.txt", sep="");
snps_location_file_name = paste(base.dir, "/2016-05-16 BXD.geno.loc.eqtl.for.lung.txt", sep="");
# Gene expression file name
expression_file_name = paste(base.dir, "/2016-05-16 mouse.lung.expression.eqtl.txt", sep="");
gene_location_file_name = paste(base.dir, "/2016-05-16 lung.gene.loc.txt", sep="");
# Covariates file name
# Set to character() for no covariates
covariates_file_name = character() ;

# Output file name
output_file_name_cis = tempfile();
output_file_name_tra = tempfile();

# Only associations significant at this level will be saved
pvOutputThreshold_cis = 1;
pvOutputThreshold_tra = 0.000000000000005;

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();
# errorCovariance = read.table("Sample_Data/errorCovariance.txt");
# Distance for local gene-SNP pairs
cisDist = 1e6;

## Load genotype data
snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;
snps$fileSkipColumns = 1;
snps$fileSliceSize = 2000;
snps$LoadFile(SNP_file_name);

## Load gene expression data
gene = SlicedData$new();
gene$fileDelimiter = "\t";
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;
gene$fileSkipColumns = 1;
gene$fileSliceSize = 2000;
gene$LoadFile(expression_file_name);

## Load covariates
cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

## Run the analysis

snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
head(genepos)

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  output_file_name = output_file_name_tra,
  pvOutputThreshold = pvOutputThreshold_tra,
  useModel = useModel,
  errorCovariance = numeric(),
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis,
  snpspos = snpspos,
  genepos = genepos,
  cisDist = cisDist,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_cis);
## Results:
cat('Analysis done in:', me$time.in.sec, ' seconds', '\n')
cat('Detected local eQTLs:','\n')
cis.eqtls<-me$cis$eqtls
head(cis.eqtls)
dim(cis.eqtls)
cis.eqtls$beta_se <-cis.eqtls$beta/cis.eqtls$statistic
write.table(cis.eqtls,file="mouselung.cis.1M.eqtls.txt", sep="\t", row.names=FALSE, quote=FALSE)


