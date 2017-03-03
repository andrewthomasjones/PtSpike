#STAT3030 Prac 2
#Soln Script
#
##load some libraries we will need

#good for reshaping data frames
#install.packages('reshape2')
library(reshape2)
#for plotting
library(lattice)
library(ggplot2)
#has special method needed in Q 2
#install.packages('EMMIXcontrasts')#only needed if not already installed
library(EMMIXcontrasts)

#Question 2
#A)
########################################################################################################
#read in data. 
goldenspike=read.csv("goldenspike.csv")

#have  alook at data
head(goldenspike)

#add extra column to goldenspike make clearer which are null genes
#new column "isNull" is TRUE when gene is a null gene
goldenspike$isNull=abs(goldenspike$fold)==1


#do t-test function
row_t.test=function(row,A_cols,B_cols){
  #var.equal=TRUE to pool variance
  test=t.test(row[A_cols], row[B_cols],var.equal=TRUE)
  return(abs(test$statistic))
}


#t test on golden spike
goldenspike$tscore=apply(goldenspike, 1, row_t.test, 1:3, 4:6)
tOrder=order(goldenspike$tscore,decreasing = TRUE)
gsSortTscore= goldenspike[tOrder,]
head(gsSortTscore,50)

#cumulative sum 
gsSortTscore$numNull= cumsum(gsSortTscore$isNull)

#plot
#xyplot(gsSortTscore$numNull[1:1000] ~ 1:1000)
qplot(x=1:1000, y=gsSortTscore$numNull[1:1000],  geom='point', xlab="N", ylab="Number of null genes", main="t-score method")

#part B
#need as matrix, only relevant columns
goldMat=as.matrix(goldenspike[,1:6])

#g is number of clusters
#debug=0 turns off verbose output
#itmax=1000,epsilon=1e-4 are stop conditions
#n1=3,n2=3 sample sizes
#ncov=3,nvcov=1 covariance structure
goldEmmix=emmixwire(goldMat,g=5,ncov=4,nvcov=1,n1=3,n2=3, debug=1,itmax=1000,epsilon=1e-4)
goldenspike$contrast=abs(scores.wire(goldEmmix))

#sort
contrastOrder=order(goldenspike$contrast,decreasing = TRUE)
gsSortContrast= goldenspike[contrastOrder,]

#get number of nulls,#cumulative sum 
gsSortContrast$numNull= cumsum(gsSortContrast$isNull)
gsSortContrast$numNull[1000]
#xyplot(gsSortContrast$numNull[1:1000] ~ 1:1000)
qplot( x=1:1000, y=gsSortContrast$numNull[1:1000],  geom='point', xlab="N", ylab="Number of null genes", main="EMMIX-Contrast method")

#plot together for niceness
combined = data.frame(n = 1:nrow(goldenspike), contrast = gsSortContrast$numNull, tscore = gsSortTscore$numNull)

combined2=melt(combined, id='n')
combined3=subset(combined2, n<=1000)

p=ggplot(data=combined3, aes(x=n, y=value, colour=variable)) + geom_point()
p=p+scale_colour_discrete(name = "Method")+ scale_x_continuous(name = "N")+ scale_y_continuous(name = "Number of null genes") 
p

#################### EMMIXwire function parameters ####################
#g is number of clusters - how many types of genes? How about 5? whats sensible for the problem?
#n1=7,n2=8 sample sizes for two sample groups 

#ncov=1,nvcov=1 covariance structure for the model- set to simplest option
#type help(emmixwire) for details - perhaps this is suboptimal

#debug=0 turns off verbose output
#itmax=1000,epsilon=1e-4 are stop conditions.
#itmax is maximum iterations. stops when change between iterations < epsilon
#############################################################


#Q2 - C
#how is the data distributed  

#### need to refomr data a bit to plot nicely

# is fold neg or pos
goldenspike$isNeg=(goldenspike$fold)<0
goldenspike$geneNum = 1:nrow(goldenspike)
goldenspike2 = melt(goldenspike, id=c("geneNum", "isNull", "isNeg", "fold", "tscore", "contrast"),variable.name ="group", value.name = "geneLevel")
#levels
levels(goldenspike2$group)
levels(goldenspike2$group)=c("G1","G1","G1","G2","G2","G2")

#tidy names on other vars too
goldenspike2$isNeg=factor(goldenspike2$isNeg)
goldenspike2$isNull=factor(goldenspike2$isNull)
levels(goldenspike2$isNeg)=c("Pos","Neg")
levels(goldenspike2$isNull)=c("NonNull","Null")
#########################
#now we are ready to plot#
p=ggplot(data=goldenspike2, aes(colour=interaction(group,isNeg,isNull), fill=interaction(group,isNeg,isNull) ))+geom_density(aes(x=geneLevel), alpha=0.4)
p+scale_colour_discrete(guide=FALSE)+scale_fill_discrete(name="Group.Fold.Null")

#additional plots may help here
# is it normal???

##################################################################
#Q2 - D
#do a function, in this case Var by groups
by(goldenspike2$geneLevel,goldenspike2$group, var)
# in this case mean by (isgroup and isNull) combined -  4 cases
by(goldenspike2$geneLevel,list(goldenspike2$group,goldenspike2$isNull),mean)

#do we want to take mean of each group first before combining

#how does the normality and variance differ among groups and null vs non null? + or - fold change?





