# Probeset normalization
#
# Qianqian Zhu, 8/2010
# 

normalize.ExpressionSet.methods = c(normalize.ExpressionSet.methods, "constant", "vsn", "loess", "group")

normalize.ExpressionSet.constant = function (eset, refindex = 1, FUN = mean, na.rm = TRUE) {
	n <- ncol(eset)
	if (!(refindex %in% 1:n))  stop("invalid reference index for normalization")
	refconstant <- FUN(exprs(eset)[, refindex], na.rm = na.rm)
	for (i in (1:n)[-refindex]) {
		m <- normalize.constant(exprs(eset)[, i], refconstant, FUN = FUN, na.rm = na.rm)
		exprs(eset)[, i] <- m
	}
	return(eset)
}

normalize.ExpressionSet.vsn = function(eset, reference, strata = NULL, subsample = nrow(eset), subset, log2scale = TRUE, log2asymp = FALSE, ...) {
        if (is.na(log2scale) || is.na(log2asymp) || (log2scale && log2asymp)) stop("'log2asymp' and 'log2scale' must not both be TRUE, and not be NA.")
        if (!missing(subset)) stop("Currently can't handle subset!") 
        vsn2res = vsn2(exprs(eset), reference = reference, returnData = FALSE, subsample = subsample, ...)
        trsfx = predict(vsn2res, newdata = exprs(eset), log2scale = log2scale)
        if (log2asymp) trsfx = (trsfx - log(2))/log(2)
        exprs(eset) <- 2^trsfx
        return(eset)
}

normalize.ExpressionSet.loess = function (eset, transfn = c("none", "log", "antilog"), span=1, ...) {
	print(span)
	transfn <- match.arg(transfn)
	if (transfn == "none") {
	    exprs(eset) <- normalize.loess(exprs(eset), span=span, ...)
	}
	else if (transfn == "antilog") {
	    exprs(eset) <- log2(normalize.loess(2^exprs(eset), span=span, ...))
	}
	else {
	    exprs(eset) <- 2^(normalize.loess(log2(exprs(eset)), span=span, ...))
	}
	return(eset)
}

normalize.ExpressionSet.group = function(eset, group, nor.method) {
	print(group)
	print(nor.method)
	samples = sampleNames(eset)
	tmp = c()
	len = length(group)
	tmp.samples = c()
	for (ngroup in 1:len) {
		tmp = cbind(tmp, exprs(normalize(eset[,group[[ngroup]]], method=nor.method)))
		tmp.samples = c(tmp.samples, samples[group[[ngroup]]])
	}
	colnames(tmp) = tmp.samples
	# make sure the sample orders are the same before and after normalization
	ind = match(samples, tmp.samples)	
	exprs(eset) = tmp[,ind]
	return(eset)
}


