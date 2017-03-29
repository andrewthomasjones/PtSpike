#install.packages(c('EMMIXcontrasts', 'ggplot2', 'reshape')) #only needed if not already installed
#load packages
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("impute")
#biocLite("samr")
#biocLite("edge")

#required
library(EMMIXcontrasts)
library(lattice)
library(edge)
library(samr)
library(limma)

#optional - for nicer plots 
library(ggplot2)
library(reshape)



#read data, needs to be in same folder
load("./platspike.Rdata")
platspike<-fit2

load("./platspike3.Rdata")
platspike<-fit3

platspike$isNull<-!platspike$truth
#add extra column to platspike make clearer which are null genes
#platspike$isNull<-abs(platspike$fold)==1


fit <- lmFit(platspike[,1:18],design=c(rep(0,9),rep(1,9)))
fit <- eBayes(fit)


gold2<-as.matrix(platspike[,1:18])
colnames(gold2)<-NULL
stud<-build_study(gold2,grp=as.factor(c(rep(0,9),rep(1,9))), sampling = "static")
t1<-odp(stud)


platspike$odpscore<-qvalueObj(t1)$stat
odpOrder<-order(platspike$odpscore,decreasing = TRUE)
odpSortTscore<- platspike[odpOrder,]

n_genes<-length(odpSortTscore$isNull)
#cumulative sum 
odpSortTscore$numNull<- cumsum(odpSortTscore$isNull)
odpSortTscore$propTrue<- (cumsum(odpSortTscore$isNull))/(1:n_genes)


#add extra column to platspike make clearer which are null genes
platspike$isNull<-abs(platspike$fold)==1


row_t.test<-function(row,A_cols,B_cols){
  #var.equal=TRUE to pool variance
  test<-t.test(row[A_cols], row[B_cols],var.equal=TRUE)
  return(abs(test$statistic))
}
platspike$limmatscore<-abs(fit$F)
gslimmaOrder<-order(platspike$limmatscore,decreasing = TRUE)
gslimmaSortTscore<- platspike[gslimmaOrder,]

n_genes<-length(gslimmaSortTscore$isNull)
#cumulative sum 
gslimmaSortTscore$numNull<- cumsum(gslimmaSortTscore$isNull)
gslimmaSortTscore$propTrue<- (cumsum(gslimmaSortTscore$isNull))/(1:n_genes)





row_t.test<-function(row,A_cols,B_cols){
  #var.equal=TRUE to pool variance
  test<-t.test(row[A_cols], row[B_cols],var.equal=TRUE)
  return(abs(test$statistic))
}

platspike$tscore<-apply(platspike, 1, row_t.test, 1:9, 10:18)
gstOrder<-order(platspike$tscore,decreasing = TRUE)
gsSortTscore<- platspike[gstOrder,]
head(gsSortTscore,50)

n_genes<-length(gsSortTscore$isNull)
#cumulative sum 
gsSortTscore$numNull<- cumsum(gsSortTscore$isNull)
gsSortTscore$propTrue<- (cumsum(gsSortTscore$isNull))/(1:n_genes)

qplot(x=1:1000, y=gsSortTscore$propTrue[1:1000],  geom=c( 'line'), xlab="N", ylab="Number of null genes", main="t-score method")


#need as matrix
goldMat<-as.matrix(platspike[,1:18])

#g is number of clusters
#debug=0 turns off verbose output
#itmax=1000,epsilon=1e-4 are stop conditions
#n1=3,n2=3 sample sizes
#ncov=3,nvcov=1 covariance structure
goldEmmix<-emmixwire(goldMat,g=4,ncov=4,nvcov=1,n1=9,n2=9, debug=1,itmax=1000,epsilon=1e-4)
platspike$contrast<-abs(scores.wire(goldEmmix))

#sort
contrastOrder<-order(platspike$contrast,decreasing = TRUE)
gsSortContrast<- platspike[contrastOrder,]

#get number of nulls,#cumulative sum 
gsSortContrast$numNull<- cumsum(gsSortContrast$isNull)
gsSortContrast$propTrue<- (cumsum(gsSortContrast$isNull))/(1:n_genes)
gsSortContrast$numNull[1000]
#xyplot(gsSortContrast$numNull[1:1000] ~ 1:1000)
qplot( x=1:1000, y=gsSortContrast$propTrue[1:1000],  geom='line', xlab="N", ylab="Number of null genes", main="EMMIX-Contrast method")

########################################################3

SAM_s0<-SAM(x=as.matrix(platspike[,1:18]), y=c(rep(1,9), rep(2,9)), resp.type="Two class unpaired", s0.perc=-1,nperms=100)
all_genes<-rbind(SAM_s0$siggenes.table$genes.up, SAM_s0$siggenes.table$genes.lo)
#sort
SAM0Order<-order(as.numeric(all_genes[,3]),decreasing = TRUE)
SAM0Sort<- platspike[as.numeric(all_genes[SAM0Order,2]),]
#get number of nulls,#cumulative sum 
SAM0Sort$numNull<- cumsum(SAM0Sort$isNull)

SAM0Sort$propTrue<- (cumsum(SAM0Sort$isNull))/(1:(nrow(SAM0Sort)))


temp<-array(0,n_genes)
temp[1:nrow(SAM0Sort)] <-SAM0Sort$propTrue

SAM0Sort_padded<-temp

#################################################################################
SAM_std<-SAM(x=as.matrix(platspike[,1:18]), y=c(rep(1,9), rep(2,9)), resp.type="Two class unpaired",nperms=100)
all_genes_2<-rbind(SAM_std$siggenes.table$genes.up, SAM_std$siggenes.table$genes.lo)
#sort
SAMOrder<-order(as.numeric(all_genes_2[,3]),decreasing = TRUE)
SAMSort<- platspike[as.numeric(all_genes_2[SAMOrder,2]),]
#get number of nulls,#cumulative sum 
SAMSort$numNull<- cumsum(SAMSort$isNull)
SAMSort$propTrue<- (cumsum(SAMSort$isNull))/(1:(nrow(SAMSort)))


temp<-array(0,n_genes)
temp[1:nrow(SAMSort)] <-SAMSort$propTrue

SAMSort_padded<-temp


#plot together for niceness
combined <- data.frame(n = 1:nrow(platspike), contrast = gsSortContrast$propTrue, tscore = gsSortTscore$propTrue, SAM_s0=SAM0Sort_padded, SAM=SAMSort_padded, limma=gslimmaSortTscore$propTrue, odp=odpSortTscore$propTrue)

combined2<-melt(combined, id='n')
combined3<-subset(combined2, n<=1000)

p<-ggplot(data=combined3, aes(x=n, y=value, colour=variable)) + geom_line()
p<-p+scale_colour_discrete(name = "Method")+ scale_x_continuous(name = "N")+ scale_y_continuous(name = "Prop True Nulls") 
p


