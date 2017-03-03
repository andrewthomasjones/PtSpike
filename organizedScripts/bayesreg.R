#   Copyright (C) 2004-2005 Suman Sundaresh, Computer Science Department
#   University of California, Irvine (suman@uci.edu)
#   Institute for Genomics and Bioinformatics
#   Parts of this code obtained or adapted from "hdarray" (Tony D. Long)

#   This code is free for academic, non-commercial, research use only. If you use this library, please
#   cite the reference (Baldi and Long, 2001) - details are on the Cyber-T
#   Help webpage. Also, please keep all headers in this code as is.
#   For commercial licenses, please contact Pierre Baldi (pfbaldi@uci.edu).


#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

#   bayesreg: v1.0beta

bayesT <- function (aData, numC, numE, ppde=TRUE, betaFit=1, bayes=TRUE, winSize=101, conf=10){
	if ((ceiling((winSize-1)/2))!=((winSize-1)/2))
		stop("ERROR: winSize must be an odd number.")

	numGene<- nrow(aData)
	#compute number of valid entries for each gene
	nC<- apply(aData[, 1:numC], 1, function(x) sum(!is.na(x)))
	nE<- apply(aData[, (numC+1):(numC+numE)], 1, function(x) sum(!is.na(x)))
	#compute means for valid entries
	meanC<- apply(aData[, 1:numC], 1, function(x) 
							if (sum(!is.na(x)))
								mean(x[!is.na(x)])
							else NA)
	meanE<- apply(aData[, (numC+1):(numC+numE)], 1, function(x) 
							if (sum(!is.na(x)))
								mean(x[!is.na(x)])
     							else NA)
	stdC<- apply(aData[, 1:numC], 1, function(x) 
							if (sum(!is.na(x)) > 1)
								sqrt(var(x[!is.na(x)]))
							else NA)
	stdE<- apply(aData[, (numC+1):(numC+numE)], 1, function(x) 
							if (sum(!is.na(x)) > 1)
								sqrt(var(x[!is.na(x)]))
							else NA)
	if (bayes){
		temp <- runavg(stdC[!is.na(stdC)][order(meanC[!is.na(stdC)])],((winSize-1)/2))	
		temp <- temp[rank(meanC[!is.na(stdC)])]	
		rasdC <- rep(NA, numGene)
		rasdC[!is.na(stdC)]<- temp

		temp <- runavg(stdE[!is.na(stdE)][order(meanE[!is.na(stdE)])],((winSize-1)/2))	
		temp <- temp[rank(meanE[!is.na(stdE)])]	
		rasdE <- rep(NA, numGene)
		rasdE[!is.na(stdE)]<- temp

		#compute bayes sd
		bayesSDC<- sqrt((conf * rasdC^2 + (nC - 1) * stdC^2)/(conf + nC - 2))
		bayesSDE<- sqrt((conf * rasdE^2 + (nE - 1) * stdE^2)/(conf + nE - 2))
		sumStats<- cbind(nC,nE,meanC,meanE,bayesSDC,bayesSDE)
		ttest<- t(apply(sumStats, 1, function(x) tstat(x)))
		colnames(ttest)=c("bayesT","bayesDF","varRatio")
		#change Bayes degree of freedom to reflect pseudo-counts
		ttest[,2]<-ttest[,2]+2*conf-2
	}
	else{
		sumStats<- cbind(nC,nE,meanC,meanE,stdC,stdE)
		ttest<- t(apply(sumStats, 1, function(x) tstat(x)))
		colnames(ttest)=c("T","DF","varRatio")
	}
	pVal <- 1 - pf(ttest[,1]^2, 1, ttest[,2])
	#calculate fold change stuff
     	fold <- - (meanC/meanE) * ((meanE/meanC) < 1) + (meanE/meanC) * ((meanC/meanE) < 1)

	#posterior probability of differential expression (ppde) test if wanted
	if(ppde){
		ppdeVal <- ppdeMix(as.matrix(pVal),betaFit)
		#if (betaFit==1)
			colnames(ppdeVal)=c("ppde<(p)","ppde(p)","ROC(x)","ROC(y)")
		#else
			#colnames(ppdeVal)=c("ppde<(p)","ppde(p)")
	}
	if (bayes){
		objBayes<- cbind(aData, nC, nE, meanC, meanE, stdC, stdE, fold, rasdC, rasdE, bayesSDC, bayesSDE,ttest, 
pVal)
		rm(nC,nE,meanC,meanE,stdC,stdE,fold,rasdC,rasdE,bayesSDC, bayesSDE, ttest,pVal,temp)
	}
	else{
		objBayes<- cbind(aData, nC, nE, meanC, meanE, stdC, stdE, fold, ttest, pVal)
		rm(nC,nE,meanC,meanE,stdC,stdE,fold,ttest,pVal)
	}
	if (ppde){
		objBayes<-cbind(objBayes, ppdeVal)
		rm(ppdeVal)
	}
	return(objBayes)
}

bayesT.pair <-  function (aData, numR, doLog=TRUE, ppde=TRUE, betaFit=1, bayes=TRUE, winSize=101, conf=10){
# NOTE: aData needs to be log-transformed
# estExp is the last column of aData and it represents
# estimated expression for that gene 
# A good value for 'estExp; would be the mean of the log of the
# "extimated expression level" over both treatments (i.e., control and
# experimentals) and replicates, where the estimated expression level for each
# gene/treatment/replicate is given as a fraction of total expression over
# all genes for that treatment/replicate.

	if ((ceiling((winSize-1)/2))!=((winSize-1)/2))
		stop("ERROR: winSize must be an odd number.")

	estExp=numR+1
	numGene<- nrow(aData)

	if (doLog==TRUE){
		cat("NOTE: The data has to be ln-transformed. The 'doLog' option is set to true,\n") 
		cat("      so performing ln-transformation on the data...\n")
		cat("      (If your data is already ln-transformed, please set 'doLog' to FALSE)\n")
		logData<-aData[,1:numR]
		logData[logData<=0]=NA
		logData<-log(logData)
		logData<-cbind(logData,aData[,estExp])
		colnames(logData)[numR+1]="estExpr"
		aData<-logData
		rm(logData)
	}
	
	#compute number of valid entries for each gene
	nR<- apply(aData[, 1:numR], 1, function(x) sum(!is.na(x)))
	#compute means for valid entries
	meanR<- apply(aData[, 1:numR], 1, function(x) 
							if (sum(!is.na(x)))
								mean(x[!is.na(x)])
							else NA)
	stdR<- apply(aData[, 1:numR], 1, function(x) 
							if (sum(!is.na(x)) > 1)
								sqrt(var(x[!is.na(x)]))
							else NA)
	if (bayes){
		index.col<- aData[,estExp]
		temp <- runavg(stdR[!is.na(stdR)][order(index.col[!is.na(stdR)])],((winSize-1)/2))	
		temp <- temp[rank(index.col[!is.na(stdR)])]	
		rasdR <- rep(NA, numGene)
		rasdR[!is.na(stdR)]<- temp

		#compute bayes sd
		bayesSD<- sqrt((conf * rasdR^2 + (nR - 1) * stdR^2)/(conf + nR - 2))
		ttest<- sqrt(nR) * (meanR/bayesSD)
		#colnames(ttest)=c("bayesT")

		#change Bayes degree of freedom to reflect pseudo-counts
		pVal <- 1 - pf(ttest^2, 1, nR + conf - 2)
	}
	else{
		ttest <- sqrt(nR) * (meanR/stdR)
		#colnames(ttest)=c("T")
		pVal <- 1 - pf(ttest^2, 1, nR - 1)
	}
	#posterior probability of differential expression (ppde) test if wanted
	if(ppde){
		ppdeVal <- ppdeMix(as.matrix(pVal),betaFit)
		colnames(ppdeVal)=c("ppde<(p)","ppde(p)","ROC(x)","ROC(y)")
	}
	if (bayes){
		objBayes<- cbind(aData, nR, meanR, stdR, rasdR, bayesSD, ttest, pVal)
		rm(nR,meanR,stdR,rasdR,bayesSD,ttest,pVal,temp)
	}
	else{
		objBayes<- cbind(aData, nR, meanR, stdR, ttest, pVal)
		rm(nR,meanR,stdR,ttest,pVal)
	}
	if (ppde){
		objBayes<-cbind(objBayes, ppdeVal)
		rm(ppdeVal)
	}
	return(objBayes)
}

tstat <- function (sumStats)
{	#number of replicates, means and std deviations
	n1=sumStats[1];  n2=sumStats[2]
	m1=sumStats[3];  m2=sumStats[4]
	sd1=sumStats[5]; sd2=sumStats[6]
	#  do the two sample t-test
	tt <- -(m1 - m2)/sqrt((((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2)/(n1 +
             n2 - 2)) * ((n1 + n2)/(n1 * n2)))
	dft <- n1 + n2 - 2
	rvar <- max((sd1^2)/(sd2^2), (sd2^2)/(sd1^2))
	as.vector(c(tt, dft, rvar))
}

ppdeMix <- function(pVal,n){
#Implements the method described in Allison et. al. 
#Computational Statistics & Data Analysis, 39:1-20 (2002) 

# The following objects are suggested objects for input into the mixture fitting
# They are lower and upper parameter bounds for the mixture parameters for a one beta,
# two beta, and three beta mixture, respectively.  The order of the parameter bounds is
# always, lambda0, lambda1, etc, then r1, s1, r2, s2, etc., depending on whether there
# are 1, 2, or 3 beta distributions.
# Note: beta parameters are capped at 170 since this value in the distribution,
# gamma(170) results in Splus reaching infinity limits.

	low1 <- c(.001,.01,.01)
	up1 <- c(1,170,170)
	low2 <- c(.001,.001,.01,.01,.01,.01)
	up2 <- c(1,1,170,170,170,170)
	low3 <- c(.001,.001,.001,.01,.01,.01,.01,.01,.01)
	up3 <- c(1,1,1,170,170,170,170,170,170)

# The following are suggested starting values for the mixture fitting
# The starting values are for parameters in the same order as above

	p01 <- c(.6,1,4)
	p02 <- c(.6,.1,1,3,1,5)
	p03 <- c(.8,.06,.1,.4,40,.4,55,1,6)

	if (ncol(pVal)==2){
		Pdata <- pVal[,2]
		Pdata <- Pdata[!is.na(Pdata)]
	}
	else
		Pdata <- pVal[!is.na(pVal)]
	Pdata<-adjPVal(Pdata)

	postP1 <- rep(NA, nrow(pVal))
	postP2 <- rep(NA, nrow(pVal))
	if(n==1){
		parms <- mle.mix(p01,mix.obj1,low1,up1,Pdata)
		L0 <- parms[1]
		L1 <- 1 - parms[1]
		r1 <- parms[2]
		s1 <- parms[3]
		cat("PPDE mixture model parameters\n")
		cat("-----------------------------\n")
		cat(paste("Lambda0 =", round(L0,5), "\n", sep=" "))
		cat(paste("Lambda1 =", round(L1,5), "\n", sep=" "))
		cat(paste("r       =", round(r1,5), "\n", sep=" "))
		cat(paste("s       =", round(s1,5), "\n", sep=" "))
   		numer1 <- L1 * pbeta(Pdata,r1,s1)
   		denom1 <- numer1 + L0 * Pdata
   		numer2 <- L1 * dbeta(Pdata,r1,s1)
   		denom2 <- numer2 + L0

		rocX <- rep(NA, nrow(pVal))
		rocY <- rep(NA, nrow(pVal))
		if (ncol(pVal)==2){
   			rocX[!is.na(pVal[,2])]<-Pdata
   			rocY[!is.na(pVal[,2])]<-pbeta(Pdata,r1,s1)
		}
		else{
   			rocX[!is.na(pVal)]<-Pdata
   			rocY[!is.na(pVal)]<-pbeta(Pdata,r1,s1)
		}
	}
	else if(n==2){
		parms <- mle.mix(p02,mix.obj2,low2,up2,Pdata)
		L0 <- parms[1]
		L1 <- parms[2]*(1-parms[1])
		L2 <-(1-parms[1])*(1-parms[2])
		r1 <- parms[3]
		s1 <- parms[4]
		r2 <- parms[5]
		s2 <- parms[6]
		cat("PPDE mixture model parameters\n")
		cat("-----------------------------\n")
		cat(paste("Lambda0 =", round(L0,5), "\n", sep=" "))
		cat(paste("Lambda1 =", round(L1,5), "\n", sep=" "))
		cat(paste("Lambda2 =", round(L2,5), "\n", sep=" "))
		cat(paste("r1      =", round(r1,5), "\n", sep=" "))
		cat(paste("s1      =", round(s1,5), "\n", sep=" "))
		cat(paste("r2      =", round(r2,5), "\n", sep=" "))
		cat(paste("s2      =", round(s2,5), "\n", sep=" "))
		numer1 <- L1 * pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2)
		denom1 <- numer1 + L0 * Pdata
		numer2 <- L1 * dbeta(Pdata,r1,s1) + L2 * dbeta(Pdata,r2,s2)
		denom2 <- numer2 + L0

		rocX <- rep(NA, nrow(pVal))
		rocY <- rep(NA, nrow(pVal))
		if (ncol(pVal)==2){
   			rocX[!is.na(pVal[,2])]<-Pdata
   			rocY[!is.na(pVal[,2])]<-(L1*pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2))/(L1+L2)
		}
		else{
   			rocX[!is.na(pVal)]<-Pdata
   			rocY[!is.na(pVal)]<-(L1*pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2))/(L1+L2)
		}
	}
	else if(n==3){
		parms <- mle.mix(p03,mix.obj3,low3,up3,Pdata)
		L0 <- parms[1]
		L1 <- parms[2]*(1-parms[1])
		L2 <- parms[3]*(1-parms[1])*(1-parms[2])
		L3 <- (1-parms[3])*(1-parms[1])*(1-parms[2])
		r1 <- parms[4]
		s1 <- parms[5]
		r2 <- parms[6]
		s2 <- parms[7]
		r3 <- parms[8]
		s3 <- parms[9]
		cat("PPDE mixture model parameters\n")
		cat("-----------------------------\n")
		cat(paste("Lambda0 =", round(L0,5), "\n", sep=" "))
		cat(paste("Lambda1 =", round(L1,5), "\n", sep=" "))
		cat(paste("Lambda2 =", round(L2,5), "\n", sep=" "))
		cat(paste("Lambda3 =", round(L3,5), "\n", sep=" "))
		cat(paste("r1      =", round(r1,5), "\n", sep=" "))
		cat(paste("s1      =", round(s1,5), "\n", sep=" "))
		cat(paste("r2      =", round(r2,5), "\n", sep=" "))
		cat(paste("s2      =", round(s2,5), "\n", sep=" "))
		cat(paste("r3      =", round(r3,5), "\n", sep=" "))
		cat(paste("s3      =", round(s3,5), "\n", sep=" "))
		numer1 <- L1 * pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2) + L3 * pbeta(Pdata,r3,s3)
		denom1 <- numer1 + L0 * Pdata
		numer2 <- L1 * dbeta(Pdata,r1,s1) + L2 * dbeta(Pdata,r2,s2) + L3 * dbeta(Pdata,r3,s3)
		denom2 <- numer2 + L0

		rocX <- rep(NA, nrow(pVal))
		rocY <- rep(NA, nrow(pVal))
		if (ncol(pVal)==2){
   			rocX[!is.na(pVal[,2])]<-Pdata
   			rocY[!is.na(pVal[,2])]<-(L1*pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2) + L3 * 
pbeta(Pdata,r3,s3))/(L1+L2+L3)
		}
		else{
   			rocX[!is.na(pVal)]<-Pdata
   			rocY[!is.na(pVal)]<-(L1*pbeta(Pdata,r1,s1) + L2 * pbeta(Pdata,r2,s2) + L3 * 
pbeta(Pdata,r3,s3))/(L1+L2+L3)
		}
	}
	else cat("Error: Must fit 1,2, or 3 beta's (change betaFit)\n")
	if (ncol(pVal)==2){
		postP1[!is.na(pVal[,2])] <- numer1/denom1
		postP2[!is.na(pVal[,2])] <- numer2/denom2
	}
	else{
		postP1[!is.na(pVal)] <- numer1/denom1
		postP2[!is.na(pVal)] <- numer2/denom2
	}
	#if (n==1)
		postP<- cbind(postP1,postP2,rocX,rocY)
	#else
	#	postP<- cbind(postP1,postP2)
	if (ncol(pVal)==2)
		postP <- cbind(pVal, postP)
	return(postP)
}

mix.obj1 <- function(p,x){
# This object is the objective function for a mixture of a uniform
# and one beta distribution
	e <- p[1] + (1-p[1])*(gamma(p[2]+p[3])/(gamma(p[2])*gamma(p[3])))*x^(p[2]-1)*(1-x)^(p[3]-1)
	sle <- -sum(log(e))
	sle
}

mix.obj2 <- function(p,x){
# This object is the objective function for a mixture of a uniform
# and two beta distributions
	e <- p[1] + (p[2]*(1-p[1]))*(gamma(p[3]+p[4])/(gamma(p[3])*gamma(p[4])))*x^(p[3]-1)*(1-x)^(p[4]-1) +
   		((1-p[1])*(1 - p[2]))*(gamma(p[5]+p[6])/(gamma(p[5])*gamma(p[6])))*x^(p[5]-1)*(1-x)^(p[6]-1)
	sle <- -sum(log(e))
	sle
}

mix.obj3 <- function(p,x){
# This object is the objective function for a mixture of a uniform
# and three beta distributions
	e <- p[1] + (p[2]*(1-p[1]))*(gamma(p[4]+p[5])/(gamma(p[4])*gamma(p[5])))*x^(p[4]-1)*(1-x)^(p[5]-1) +
		(p[3]*(1-p[1])*(1 - p[2]))*(gamma(p[6]+p[7])/(gamma(p[6])*gamma(p[7])))*x^(p[6]-1)*(1-x)^(p[7]-1) +
		((1-p[3])*(1-p[1])*(1 - p[2]))*(gamma(p[8]+p[9])/(gamma(p[8])*gamma(p[9])))*x^(p[8]-1)*(1-x)^(p[9]-1)
		sle <- -sum(log(e))
		sle
}

mle.mix <- function(init,mix.obj,low,up,Pdata){
# This object computes the MLE's for the mixture distribution
# Inputs are:
# init = starting values for the computations (see p01, p02, p03 above)
# mix.obj = is the objective function (i.e., see mix.obj1 etc from above)
# low and up are the parameter limits (see low1, up1 etc from above)
# Pdata = the p-values.  This should be a single vector of length
# equal to the number of genes in the study.
# Note: all inputs should be consistent with regards to the number of components
# being fitted.  That is, p01 goes with low1, up1, and mix.obj1.
# Output:  An object called mix that is the result of the fitting algorithm.
# mix$par are the MLE's of parameters in the same order as above objects.
# Note: the weighting parameters on mixture components need to be reparametrized
# to make sense.  
	mix <- optim(init,mix.obj,lower=low, upper=up, method="L-BFGS-B", x=Pdata)           
	return(mix$par)
}

runavg <- function(x,k=1) {
# x is the input vector
# n is the length of the vector
# k is the size of the running average window, the window 
#   includes the data point at that position plus k points 
#   on either side.  For example k = 3 would be a sliding
#   window of size {3+1+3}= 7 centered on each point
# r is the output vector
     	if (k <= 0)
        	 stop("Error: Window size (k) must be greater than zero")
	x<-as.array(x);
     	n <- length(x)

	r<-0
	t<-0

	for (j in 0:k) {
		t<-t+x[j+1]
	}

	l=0
	u=k

	for (i in 0:(n-1)) {
	# the use of u-l-1 instead of k handles end effects 
		r[i+1]= (t/(u-l+1))
	# update the current running total for the next position 
		if (i > k-1) {
			t=t-x[l+1]
			l=l+1
		}
		if (u < n-1) {
			u=u+1			
			t=t+x[u+1]
		}
	}

	return(as.array(r))
}

runavgPool <- function(x,k=1,start,end) {
# x is the input vector
# start and end mark the column numbers of expt or control data
# n is the length of the vector
# k is the size of the running average window, the window 
#   includes the data point at that position plus k points 
#   on either side.  For example k = 3 would be a sliding
#   window of size {3+1+3}= 7 centered on each point
# r is the output vector
     	if (k <= 0)
        	 stop("Error: Window size (k) must be greater than zero")

	x<-as.matrix(x[,start:end]);
     	n <- dim(x)[1]

	r<-0
	t<-0

	side <- (k-1)/2
	
	for (j in 1:side) {
		flatArr <- NULL
		for(i in 1:(j+side))
			flatArr <- c(flatArr,x[i,])

		flatArr<-as.array(flatArr)
		
		r[j]<-sd(flatArr)		
	}
	for (j in (side+1):(n-(side+1))) {
		flatArr <- NULL
		for (i in (j-side):(j+side)) 
			flatArr <- c(flatArr,x[i,])
		flatArr<-as.array(flatArr)	
		r[j]<-sd(flatArr)		

	}
	for (j in (n-side):n) {
		flatArr <- NULL
		for(i in (j-side):n)
			flatArr <- c(flatArr,x[i,])
		flatArr<-as.array(flatArr)	
		r[j]<-sd(flatArr)		
	}
	return(as.array(r))
}

secMin<-function(a){
#returns the next minimum after the overall min
	i=0
	a<-sort(a)
	secondMin=a[1]
	while (secondMin==a[1]){
		i=i+1
		secondMin=a[i]
	}
	return(secondMin)
}

secMax<-function(a){
#returns the next maximum after the overall max
	i=0
	a<-sort(a,decreasing=TRUE)
	secondMax=a[1]
	while (secondMax==a[1]){
		i=i+1
		secondMax=a[i]
	}
	return(secondMax)
}


adjPVal<-function(a){
#makes sure that p-values are not zero but some very small value (factor of 10 smaller than second min)
#also if p-values are 1, then they are made a very small number less than 1 but more than all others
	if(min(a)==0){
		adj<-secMin(a)/10
		a[a==0]=adj
	}
	if(max(a)==1){
		adj<-secMax(a) + ((1-secMax(a))*9/10)
		a[a==1]=adj
	}
	return(a)
}



