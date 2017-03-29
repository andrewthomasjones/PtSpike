#
# DEGAUC.R
#
# calculate partial AUC values
#
# citation: Zhu Q, Miecznikowski J, Halfon M (2010). Preferred analysis methods for Affymetrix GeneChips. II. An expanded, balanced, wholly-defined spike-in dataset. BMC Bioinformatics. (2010) 11:285.
#
# Qianqian Zhu, 8/2010


library(affy)
library(ROCR)

# Note: In ROCR, f(x) >= cutoff => class "+"; f(x) < cutoff => class "-". 
# Use abs(statistic) to evaluate DEGs: the higher value, the more likely to be DEGs.
# abs(statistic) >= cutoff => class "+" (DEGs) 
fcinfo = read.table("./Affyprobeset.fcvalue", header=F)
truth = rep(NA, nrow(fcinfo))
# empty probesets are considered as not differentially expressed probesets 
truth[fcinfo$V2!=0 & fcinfo$V2 != "MC" & fcinfo$V2!="MF" & fcinfo$V2!="NS" & fcinfo$V2!=1] = 1 # differentially expressed probesets
truth[fcinfo$V2==1] = 0 # not differentially expressed probesets
truth[fcinfo$V2==0] = 0 # not differentially expressed probesets
names(truth) = fcinfo$V1
truth = truth[!is.na(truth)]


methods = c("gcrma-rs.scaling.3.pmonly.medianpolish.none.1.cyberT", "gcrma-rs.scaling.3.pmonly.medianpolish.none.1.samr",
            "gcrma-rs.constant.3.pmonly.medianpolish.none.1.cyberT", "gcrma-rs.scaling.3.pmonly.medianpolish.none.1.limma", "gcrma-rs.constant.3.pmonly.medianpolish.none.1.limma") # 9 of the top 10 routes. Achemy method can't not be performed using this script.

#1,3 only work
for (m in methods) {
	cat(m, sep="\n")
	file = paste("test/pt", m, "RData", sep=".")
	cat(file, sep="\n")
	load(file)
	fit2<-as.data.frame(fit[,1:18])
	
	# the probesets that have expected fold change information
	index = which(is.na(match(names(fit[,1]), names(truth))))
	fit2 = fit2[-index,]
	truth.tmp = truth[match(rownames(fit2), names(truth))]
	
	fit2$truth<-truth.tmp
}
