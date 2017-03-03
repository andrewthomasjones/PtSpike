# Probe normalization using only pm probes 
#
# Qianqian Zhu, 8/2010
#
# Note:
# constantpm requires package "affy"
# scalingpm requires package "affyPLM"
# vsnpm requires package "vsn"

normalize.AffyBatch.methods = c(normalize.AffyBatch.methods, "constantpm", "scalingpm", "invariantsetpm", "loesspm", "quantilespm", "vsnpm")

normalize.AffyBatch.constantpm <- function(abatch, refindex=1, FUN=mean, na.rm=TRUE) {
	n <- length(abatch)
	Index <- unlist(indexProbes(abatch, which="pm"))
	print("Got PM indices...")
	if (! (refindex %in% 1:n)) stop("invalid reference index for normalization")
	refconstant <- FUN(intensity(abatch)[Index, refindex], na.rm=na.rm)
	normhisto <- vector("list", length = n)
	for (i in (1:n)[-refindex]) {
		m <- normalize.constant(intensity(abatch)[Index, i], refconstant, FUN=FUN, na.rm=na.rm)
		intensity(abatch)[Index, i] <- m
		myhistory <- list(name="normalized by constant", constant=attr(m,"constant"))
		attr(m,"constant") <- NULL
		normhisto[[i]] <- myhistory
	}
	attr(abatch, "normalization") <- normhisto
	return(abatch)
}

normalize.AffyBatch.scalingpm <- function (abatch, type = "pmonly",  trim = 0.02, baseline = -1, log.scalefactors = FALSE) {
	return(normalize.AffyBatch.scaling(abatch, type=type, trim=trim, baseline=baseline, log.scalefactors=log.scalefactors))
}

normalize.AffyBatch.invariantsetpm <- function (abatch, prd.td = c(0.003, 0.007), verbose = FALSE, baseline.type = c("mean", "median", "pseudo-mean", "pseudo-median"), type = "pmonly") {
	return(normalize.AffyBatch.invariantset(abatch, prd.td=prd.td, verbose=verbose, baseline.type=baseline.type, type=type))
}

normalize.AffyBatch.loesspm <- function (abatch, type = "pmonly") {
	return(normalize.AffyBatch.loess(abatch, type=type))
}

normalize.AffyBatch.quantilespm <- function (abatch, type = "pmonly") {
	return(normalize.AffyBatch.quantiles(abatch, type=type))
}

normalize.AffyBatch.vsnpm <- function(abatch, reference, strata = NULL, subsample = 30000L, subset, log2scale = TRUE, log2asymp = FALSE, ...) {
	if (is.na(log2scale) || is.na(log2asymp) || (log2scale && log2asymp)) stop("'log2asymp' and 'log2scale' must not both be TRUE, and not be NA.")
	ind = indexProbes(abatch, "pm")
	if (!missing(subset)) ind = ind[subset]
	ind = unlist(ind)
	vsn2res = vsn2(intensity(abatch)[ind, ], reference = reference, returnData = FALSE, subsample = subsample, ...)
	description(abatch)@preprocessing = c(description(abatch)@preprocessing, list(vsnReference = vsn2res))
	trsfx = predict(vsn2res, newdata = intensity(abatch)[ind,], log2scale = log2scale)
	if (log2asymp) trsfx = (trsfx - log(2))/log(2)
	intensity(abatch)[ind,] <- 2^trsfx # QQ
	return(abatch)
}
	
