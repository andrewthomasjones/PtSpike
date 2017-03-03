################### STAT3003 Semester 1, 2016 ###################
#Practical 1
#
#
#
################### LOADING DATA AND R-PACKAGES ###########################

#install package EMMIXcontrasts
install.packages('EMMIXcontrasts') #only needed if not already installed
library(EMMIXcontrasts)

#read data, needs to be in working directory
load("./hedenfalk.Rdata")
load("./goldenspike.Rdata")

#add extra column to goldenspike make clearer which are null genes
#new column "isNull" is TRUE when gene is a null gene
goldenspike$isNull<-abs(goldenspike$fold)==1

################### USEFUL FUNCTIONS ############################

#returns the pooled t-test test stat for test on Group1 and Group2
t.test(group1, group2, var.equal=TRUE)$statistic

#takes a vector or list and returns it sorted from biggest to smallest
ordering<-order(unordered, decreasing = TRUE)

#returns a dataframe reordered (by row) according to ordering
dataframe[ordering,]

# returns the cumulative sum of a vector i.e. cumsum(1:3) = 1 3 6
cumsum(x)

###################### USING EMMIXcontrasts ######################

#need data as matrix, only want first 15 columns which contain the data
hfMat<-as.matrix(hedenfalk[,1:15])

hfEmmix<-emmixwire(hfMat,g=5,ncov=1,nvcov=1,n1=7,n2=8, debug=1, itmax=1000,epsilon=1e-4)

#return absolute values of score from fitted emmix model
hfcontrast<-abs(scores.wire(hfEmmix))

#################### EMMIXwire function parameters ####################

#g is number of clusters - how many types of genes? How about 5? whats sensible for the problem?
#n1=7,n2=8 sample sizes for two sample groups 

#ncov=1,nvcov=1 covariance structure for the model- set to simplest option
#type help(emmixwire) for details - perhaps this is suboptimal

#debug=0 turns off verbose output
#itmax=1000,epsilon=1e-4 are stop conditions.
#itmax is maximum iterations. stops when change between iterations < epsilon

#############################################################