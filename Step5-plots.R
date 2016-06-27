liver.mouse.eQTL.bayesian.tau <- read.table("liver.mouse.eQTL.bayesian.tau.txt")

Fcomb = function(ps) #chi-square (Fisher, 1932, Lancaster, 1961)
{
        k = length(ps)
        temp = -2*sum(log(ps))
        pchisq(temp, 2*k, lower.tail = F)
}

Ncomb = function(ps) #normal (Liptak, 1958, Stouffer 1949)
{
        k = length(ps)
        z = qnorm((1-ps))
        Ts = sum(z)/sqrt(k) # sum(1-Phi^-1(1-p))/sqrt(k)
        pnorm(Ts, lower.tail = F)  #Same as 1-Phi

}

#META
metapval = apply(cbind(liver.mouse.eQTL.bayesian.tau$lung_pvalue, liver.mouse.eQTL.bayesian.tau$liver_pvalue), 1, Ncomb)
liver.mouse.eQTL.bayesian.tau$metapval = metapval

#Multiple posterior prob by 2 
##CHANGE?
#liver.mouse.eQTL.bayesian.tau$p.below.0 = 2*liver.mouse.eQTL.bayesian.tau$p.below.0

#MT Method
mtresults<-read.table(paste0("MTeQTLs_ASE_3c_",sebsetn,"s.txt"), header = TRUE)

minmtresults<-sapply(liver.mouse.eQTL.bayesian.tau$ensembl_id, function(x) max(mtresults[mtresults$ensembl_id == as.character(x),"marginalP.liver"]))
newmtresults = data.frame(liver.mouse.eQTL.bayesian.tau$ensembl_id, minmtresults, liver.mouse.eQTL.bayesian.tau$eqtl)
colnames(newmtresults) = c("ensembl_id", "marginalp", "eqtl")

newresults = liver.mouse.eQTL.bayesian.tau[,c("ensembl_id", "lung_pvalue", "liver_pvalue", "metapval", "p.below.0", "eqtl")]
pvals = sort(newresults[,"lung_pvalue"])
nvals = length(newresults[,"lung_pvalue"])
result.lung = result.liver = result.meta = result.pp = result.mt = matrix(0, nvals, 3)
totaltrue = sum(newresults[,"eqtl"])
totalfalse = sum(newresults[,"eqtl"]==0)
j = 1
for (i in pvals)
{
   result.lung[j,1] = i
   result.lung[j,2] = sum( newresults[newresults[,"lung_pvalue"]<i,"eqtl"])/totaltrue  #sens
   result.lung[j,3] = sum( newresults[newresults[,"lung_pvalue"]>=i,"eqtl"]==0)/totalfalse #spec

   result.liver[j,1] = i
   result.liver[j,2] = sum( newresults[newresults[,"liver_pvalue"]<i,"eqtl"])/totaltrue  #sens
   result.liver[j,3] = sum( newresults[newresults[,"liver_pvalue"]>=i,"eqtl"]==0)/totalfalse #spec

   result.mt[j,1] = i
   result.mt[j,2] = sum( newmtresults[newmtresults[,"marginalp"]<i,"eqtl"])/totaltrue  #sens
   result.mt[j,3] = sum( newmtresults[newmtresults[,"marginalp"]>=i,"eqtl"]==0)/totalfalse #spec

   result.meta[j,1] = i
   result.meta[j,2] = sum( newresults[newresults[,"metapval"]<i,"eqtl"], na.rm  = TRUE)/sum(newresults[,"eqtl"]==1, na.rm = TRUE) #sens
   result.meta[j,3] = sum( newresults[newresults[,"metapval"]>=i,"eqtl"]==0, na.rm  = TRUE)/sum(newresults[,"eqtl"]==0, na.rm = TRUE) #spec

   result.pp[j,1] = i
   result.pp[j,2] = sum( newresults[newresults[,"p.below.0"]<i,"eqtl"], na.rm  = TRUE)/sum(newresults[,"eqtl"]==1, na.rm = TRUE)  #sens
   result.pp[j,3] = sum( newresults[newresults[,"p.below.0"]>=i,"eqtl"]==0, na.rm  = TRUE)/sum(newresults[,"eqtl"]==0, na.rm = TRUE) #spec

   j = j+1
}


plot(1-result.liver[,3], result.liver[,2], xlim = c(0,1), ylim = c(0,1), xlab= "1-Specity", ylab = "Sensitivy", pch = ".")
points(1-result.pp[,3], result.pp[,2], col = "red", pch = ".")
points(1-result.lung[,3], result.lung[,2], col = "purple", pch = ".")
points(1-result.meta[,3], result.meta[,2], col = "blue", pch = ".")
points(1-result.mt[,3], result.mt[,2], col = "green", pch = ".")
legend(.7, .55, legend = c("orig - liver", "bayesian", "lung", "meta", "mt"), col = c("black", "red", "purple", "blue", "green"), pch = 1)
title(paste0("Z:abs_lung_beta; subsample:", sebsetn))


 dev.copy(pdf,"comparison.pdf")
dev.off()

library(flux)
orig_auc <- auc(1-result.liver[,3], result.liver[,2])
bayesian_auc <- auc(1-result.pp[,3], result.pp[,2])
lung_auc <- auc(1-result.lung[,3], result.lung[,2])
meta_auc <- auc(1-result.meta[,3], result.meta[,2])
mt_auc <- auc(1-result.mt[,3], result.mt[,2])
auc <- rbind(orig_auc, bayesian_auc,lung_auc,meta_auc,mt_auc)
auc <- cbind(auc, auc[ , 1]/auc["orig_auc", 1])
colnames(auc) <- c("auc", "FC")
auc
