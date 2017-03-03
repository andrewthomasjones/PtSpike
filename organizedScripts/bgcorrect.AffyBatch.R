# Background correction using gcrma
#
# Qianqian Zhu, 8/2010
#

bgcorrect.methods = c(bgcorrect.methods, "gcrma")
bg.correct.gcrma = function(abatch, ...) {
	abatch = bg.adjust.gcrma(abatch, ...)
	return(abatch)
} 
