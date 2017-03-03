# DEG testing
# 
# Qianqian Zhu, 8/2010
#
# Note: the intensities are in natural scale (not log scale)
#	and the arrays in the same condition are next to each other in eset, eg: C1, C2, C3, E1, E2, E3.
#	C, E are the two conditions.
#	bayesreg.R is required for DEG testing using cyberT.

library(limma)
library(samr)

DEGtest.cyberT = function(eset, C.ind, S.ind, ppde=F, ...) {
	all.exprs = log2(exprs(eset))
	fit = bayesT(all.exprs, numC=length(C.ind), numE=length(S.ind), ppde=ppde, ...)
	return(fit)
}

DEGtest.limma = function(eset, C.ind, S.ind) {
	all.exprs = exprs(eset)
	ind = which(all.exprs==0) # if there are some probesets having intensities 0, log2(0) = -Inf and lead to error when performing limma
	if (length(ind)>0) {
		cat("0 intensity!", sep="\n")
		all.exprs[ind] = all.exprs[ind] + 10^(-10) 
	}
	exprs(eset) = log2(all.exprs)
	design = data.frame(C=rep(1, ncol(all.exprs)), S=c(rep(0, length(C.ind)), rep(1, length(S.ind))), row.names=sampleNames(eset))
	fit = lmFit(eset, design)
	fit = eBayes(fit)
	return(fit)
}

DEGtest.samr = function(eset, C.ind, S.ind) {
	all.exprs = log2(exprs(eset))
	samples = c(rep("1", length(C.ind)), rep("2", length(S.ind)))
	data = list(x=all.exprs, y=samples, logged2=T)
	fit = samr(data, resp.type="Two class unpaired", nperms=1) 
	return(fit)
}

DEGtest.fc = function(eset, C.ind, S.ind) {
	all.exprs = exprs(eset)
	fit = rowMeans(all.exprs[,S.ind])/rowMeans(all.exprs[,C.ind])
	names(fit) = rownames(all.exprs)
	return(fit)
}
		
