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
fcinfo = read.table("Affyprobeset.fcvalue", header=F)
truth = rep(NA, nrow(fcinfo))
# empty probesets are considered as not differentially expressed probesets 
truth[fcinfo$V2!=0 & fcinfo$V2 != "MC" & fcinfo$V2!="MF" & fcinfo$V2!="NS" & fcinfo$V2!=1] = 1 # differentially expressed probesets
truth[fcinfo$V2==1] = 0 # not differentially expressed probesets
truth[fcinfo$V2==0] = 0 # not differentially expressed probesets
names(truth) = fcinfo$V1
truth = truth[!is.na(truth)]

methods = c("gcrma-rs.scaling.3.pmonly.medianpolish.vsn.1.cyberT", "none.vsn.3.pmonly.medianpolish.constant.3.samr", "gcrma-rs.scaling.3.pmonly.medianpolish.vsn.1.samr",
		"gcrma-rs.constant.3.pmonly.medianpolish.vsn.1.cyberT", "gcrma-rs.scaling.3.pmonly.medianpolish.vsn.1.limma", "gcrma-rs.constant.3.pmonly.medianpolish.vsn.1.limma",
		"none.vsn.3.pmonly.medianpolish.constant.3.limma", "none.vsn.3.pmonly.medianpolish.quantiles.1.fc", "rma.vsn.1.pmonly.medianpolish.scaling.3.samr") # 9 of the top 10 routes


fpr.stop = 0.05 # AUC is the area under the curve when fpr <= fpr.stop
write.table(t(c("n", paste("AUC", fpr.stop, sep="_"))), file="test/DEGAUC.diag", append=F, row.names=F, col.names=F, quote=F, sep="\t")
for (m in methods) {
	cat(m, sep="\n")
	file = paste("test/pt", m, "RData", sep=".")
	cat(file, sep="\n")
	load(file)
	DEGm = strsplit(m, "\\.")[[1]][8]
	if (DEGm == "cyberT") {
		statistic = fit[,which(colnames(fit)=="bayesT")]
		names(statistic) = rownames(fit)
	} else {
		if (DEGm == "limma") {
			statistic = fit$t[,2]
			names(statistic) = rownames(fit$t)
		} else {
			if (DEGm == "fc") {
				statistic = log2(fit)
				names(statistic) = names(fit)
			} else {
				if (DEGm == "samr") {
					statistic = fit$tt
				}
			}
		}
	}
	# the probesets that have expected fold change information
	index = which(is.na(match(names(statistic), names(truth))))
	statistic = statistic[-index]
	truth.tmp = truth[match(names(statistic), names(truth))]
	pred = prediction(predictions=abs(statistic), labels=truth.tmp)
	diag.stat = c(length(truth.tmp), mapply(function(x){attr(performance(pred, measure="auc", fpr.stop=x), "y.values")[[1]]}, fpr.stop))
	write.table(t(c(m, diag.stat)), file="test/DEGAUC.diag", append=T, row.names=F, col.names=F, quote=F, sep="\t")
}
	
warnings()
