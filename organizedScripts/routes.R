# routes.R
# 
# from background correction to DEG testing 
#
# citation: Zhu Q, Miecznikowski J, Halfon M (2010). Preferred analysis methods for Affymetrix GeneChips. II. An expanded, balanced, wholly-defined spike-in dataset. BMC Bioinformatics. (2010) 11:285.
# 
# Qianqian Zhu, 8/2010
 
library(affy)
library(affyPLM) # for scaling normalization, which is the same as "constant" but using trim mean; by default, type="together", trim=0.02 in function normalize.AffyBatch.scaling
normalize.AffyBatch.methods = c(normalize.AffyBatch.methods, "scaling")
library(vsn) # for vsn normalization
source("./organizedScripts/bgcorrect.AffyBatch.R") # bg.correct.gcrma  
source("./organizedScripts/normalizegroup.AffyBatch.R") # normalize.AffyBatch.group
source("./organizedScripts/normalizepm.R") # normalize between pm probes
source("./organizedScripts/expresso.R")
source("./organizedScripts/probesetnormalize.R") # probeset normalization
source("./organizedScripts/fullroute.R")
source("./organizedScripts/bayesreg.R") # cyberT
library(limma)
library(samr)
source("./organizedScripts/DEGtests.R")

methods = c("gcrma-rs.scaling.3.pmonly.medianpolish.vsn.1.cyberT", "none.vsn.3.pmonly.medianpolish.constant.3.samr", "gcrma-rs.scaling.3.pmonly.medianpolish.vsn.1.samr",
		"gcrma-rs.constant.3.pmonly.medianpolish.vsn.1.cyberT", "gcrma-rs.scaling.3.pmonly.medianpolish.vsn.1.limma", "gcrma-rs.constant.3.pmonly.medianpolish.vsn.1.limma",
		"none.vsn.3.pmonly.medianpolish.constant.3.limma", "none.vsn.3.pmonly.medianpolish.quantiles.1.fc", "rma.vsn.1.pmonly.medianpolish.scaling.3.samr") # 9 of the top 10 routes. Achemy method can't not be performed using this script.
# Note: 
# "gcrma-rs" corresponds to "gcrma-reb" in publication.
# "1" means technical normalization, "2" means conditional normalization and "3" means normalization using all arrays.

nmethods = length(methods)
print(methods)

load("test/Affydata.raw.RData")
conditions = c("A", "B")
sample.names = sampleNames(affydata)
print(sample.names)
ind = mapply(function(x){grep(x, sample.names)}, conditions)
print(ind)
affydata = affydata[,as.vector(ind)] # organize the arrays to make sure the arrays in same condition are next to each other
ind = mapply(function(x){grep(x, sampleNames(affydata))}, conditions)
mapply(fullroute, method=methods, MoreArgs=list(affydata=affydata, label="test/pt", C.ind=ind[,1], S.ind=ind[,2])) 

warnings()

