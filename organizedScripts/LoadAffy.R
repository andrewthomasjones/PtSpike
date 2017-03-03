#
# LoadAffy.R
#
# load the Platinum Spike dataset and save as a R object
#
# citation: Zhu Q, Miecznikowski J, Halfon M (2010). Preferred analysis methods for Affymetrix GeneChips. II. An expanded, balanced, wholly-defined spike-in dataset. BMC Bioinformatics. (2010) 11:285.
#
# Qianqian Zhu, 8/2010
#

library(affy)

# load CEL files
fnames = list.files(path="../cel_files", pattern="CEL", full.names=T)
nsamples = length(fnames)
affydata = ReadAffy(filenames=fnames)

# name samples
sample.name = sampleNames(affydata)

getlbl = function(i) {
        strsplit(strsplit(sample.name[i], "_")[[1]][2], "\\.")[[1]][1]
}
labels = mapply(getlbl, 1:nsamples)

# reorder samples to F1A1, F1A2, F1A3, F2A1, F2A2, F2A3, F3A1, F3A2, F3A3, F1B1, F1B2, F1B3, F2B1, F2B2, F2B3, F3B1, F3B2, F3B3
neworder = order(substr(labels, 3, 3))

labels = labels[neworder]
affydata = affydata[,neworder]

sampleNames(affydata) = labels

save(file="../test/Affydata.raw.RData", affydata, compress=T)

warnings()

