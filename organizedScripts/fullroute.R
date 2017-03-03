# From background correction to DEG testing
#
# Qianqian Zhu, 8/2010

fullroute = function(method, affydata, label, C.ind, S.ind) {
        cat(method, sep="\n")
        mets = strsplit(method, "\\.")[[1]]
        # background correction
        bg.correct = T
        bgcorrect.method = mets[1]
        bgcorrect.param = list()
        if (bgcorrect.method == "none") bg.correct = F
        if (bgcorrect.method == "gcrma-rs") {
                bgcorrect.method = "rma"
                bgcorrect.param = list(fast=F)
        }
        if (bgcorrect.method == "gcrma-lf") {
                bgcorrect.method = "rma"
                bgcorrect.param = list(affinity.source="local")
        }        
	if (bgcorrect.method == "gcrma-ls") {
                bgcorrect.method = "rma"
                bgcorrect.param = list(affinity.source="local", fast=F)
        }
        # probe normalization
        samples = sampleNames(affydata)
        nor.group = c()
        nor.group[1] = list(unique(substr(samples, 1, 3))) # only normalize between the technical replicates
        nor.group[2] = list(unique(substr(samples, 3, 3))) # normalize between the same condition: A_B
        nor.group[3] = list(unique(substr(samples, 1, 1))) # normalize using all samples
        normalize = F
        normalize.param = list()
	normalize.group = "" 
        if (mets[2] != "none") {
                normalize = T
                normalize.method = mets[2]
                normalize.group = as.numeric(mets[3])
                group = c()
                for (ng in 1:length(nor.group[[normalize.group]])) {
                        group[ng] = list(grep(nor.group[[normalize.group]][ng], samples))
                }
		# Note: Can't use method=normalize.tmp[1], otherwise in expresso function, normalize.method is overwriten by method in afbatch <- do.call("normalize", c(alist(afbatch, normalize.method), normalize.param))
		normalize.param = list(group=group, nor.method=mets[2])
                if (length(grep("rma", bgcorrect.method)) > 0) normalize.param = list(group=group, nor.method=paste(mets[2], "pm", sep="")) 
        }
        # pm correction
        if (mets[4] == "none") {
                 pmcorrect.method = NULL
        } else {
                pmcorrect.method = mets[4]
        }
        # summarization
        if (mets[5] == "none") {
                summary.method = NULL
        } else {
                summary.method = mets[5]
        }

        #eset = expresso(affydata, bg.correct=bg.correct, bgcorrect.method=bgcorrect.method, bgcorrect.param=bgcorrect.param, normalize=normalize, normalize.method=normalize.method, normalize.param=normalize.param, pmcorrect.method=pmcorrect.method, summary.method=summary.method)
        eset = expresso(affydata, bg.correct=bg.correct, bgcorrect.method=bgcorrect.method, bgcorrect.param=bgcorrect.param, normalize=normalize, normalize.method=normalize.method, pmcorrect.method=pmcorrect.method, summary.method=summary.method)
        
        if (!is.null(summary.method)) {
                if (summary.method == "medianpolish" | summary.method == "farms" | summary.method == "dfw") exprs(eset) = 2^exprs(eset) # Note: the intensities returned by medianpolish and farms are in log2 scale, so I change it back to natural scale
        }
       
        # probeset normalization
        if (mets[6] != "none") {
                if (normalize.group != as.numeric(mets[7])) {
                        normalize.group = as.numeric(mets[7])
                        group = c()
                        for (ng in 1:length(nor.group[[normalize.group]])) {
                                group[ng] = list(grep(nor.group[[normalize.group]][ng], samples))
                        }
                }
                normalize.param = list(group=group, nor.method=mets[6])
                eset = do.call("normalize", c(alist(eset, "group"), normalize.param))
        }

	output = paste(label, method, "RData", sep=".")
	# DEG test
	if (!is.na(mets[8])) {
		fit = do.call(paste("DEGtest", mets[8], sep="."), alist(eset, C.ind=C.ind, S.ind=S.ind))
		save(eset, file=paste(label, sub(paste("\\.", mets[8], "$", sep=""), "", method), "RData", sep="."), compress=T)
		save(fit, file=output, compress=T)
	} else {
        	save(eset, file=output, compress=T)
	}	
	cat(output, sep="\n")	
}

