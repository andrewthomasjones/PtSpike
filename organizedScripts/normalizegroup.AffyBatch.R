# Probe normalization using different groups of arrays
#
# Qianqian Zhu, 8/2010
# 

normalize.AffyBatch.methods = c(normalize.AffyBatch.methods, "group")

normalize.AffyBatch.group = function(abatch, group, nor.method) {
	print(group)
	print(nor.method)
	samples = sampleNames(abatch)
	tmp = c()
	len = length(group)
	for (ngroup in 1:len) {
		tmp[ngroup] = list(normalize(abatch[,group[[ngroup]]], method=nor.method))
	}
	abatch = tmp[[1]]
	if (len >= 2) {
		for (ngroup in 2:len) {
			abatch = merge(abatch, tmp[[ngroup]])
		}
	}
	# make sure the sample orders are the same before and after normalization
	return(abatch[,match(samples, sampleNames(abatch))])
}

