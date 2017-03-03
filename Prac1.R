#install.packages(c('EMMIXcontrasts', 'ggplot2', 'reshape')) #only needed if not already installed
#load packages

#required
library(EMMIXcontrasts)
library(lattice)

#optional - for nicer plots 
library(ggplot2)
library(reshape)

#read data, needs to be in same folder
load("./hedenfalk.Rdata")
load("./goldenspike.Rdata")

#add extra column to goldenspike make clearer which are null genes
goldenspike$isNull<-abs(goldenspike$fold)==1

#question 1
#part A
#do test function
row_t.test<-function(row,A_cols,B_cols){
  #var.equal=TRUE to pool variance
  test<-t.test(row[A_cols], row[B_cols],var.equal=TRUE)
  return(abs(test$statistic))
}

hedenfalk$tscore<-apply(hedenfalk, 1, row_t.test, 1:7, 8:15)
hftOrder<-order(hedenfalk$tscore,decreasing = TRUE)
hfSortTscore<- hedenfalk[hftOrder,]
head(hfSortTscore,50)

#part B
goldenspike$tscore<-apply(goldenspike, 1, row_t.test, 1:3, 4:6)
gstOrder<-order(goldenspike$tscore,decreasing = TRUE)
gsSortTscore<- goldenspike[gstOrder,]
head(gsSortTscore,50)

#Question 2

#part A

#cumulative sum 
gsSortTscore$numNull<- cumsum(gsSortTscore$isNull)
#plot
#xyplot(gsSortTscore$numNull[1:1000] ~ 1:1000)
qplot(x=1:1000, y=gsSortTscore$numNull[1:1000],  geom='point', xlab="N", ylab="Number of null genes", main="t-score method")

#part B
#need as matrix
goldMat<-as.matrix(goldenspike[,1:6])

#g is number of clusters
#debug=0 turns off verbose output
#itmax=1000,epsilon=1e-4 are stop conditions
#n1=3,n2=3 sample sizes
#ncov=3,nvcov=1 covariance structure
goldEmmix<-emmixwire(goldMat,g=4,ncov=4,nvcov=1,n1=3,n2=3, debug=1,itmax=1000,epsilon=1e-4)
goldenspike$contrast<-abs(scores.wire(goldEmmix))

#sort
contrastOrder<-order(goldenspike$contrast,decreasing = TRUE)
gsSortContrast<- goldenspike[contrastOrder,]

#get number of nulls,#cumulative sum 
gsSortContrast$numNull<- cumsum(gsSortContrast$isNull)
gsSortContrast$numNull[1000]
#xyplot(gsSortContrast$numNull[1:1000] ~ 1:1000)
qplot( x=1:1000, y=gsSortContrast$numNull[1:1000],  geom='point', xlab="N", ylab="Number of null genes", main="EMMIX-Contrast method")

#plot together for niceness
combined <- data.frame(n = 1:nrow(goldenspike), contrast = gsSortContrast$numNull, tscore = gsSortTscore$numNull)

combined2<-melt(combined, id='n')
combined3<-subset(combined2, n<=1000)

p<-ggplot(data=combined3, aes(x=n, y=value, colour=variable)) + geom_point()
p<-p+scale_colour_discrete(name = "Method")+ scale_x_continuous(name = "N")+ scale_y_continuous(name = "Number of null genes") 
p



