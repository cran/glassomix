\name{summary.glasso.mix}
\alias{summary.glasso.mix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summary function according to function ``glasso.mix''
}
\description{
Reduced summary of the result according to ``glasso.mix''
}
\usage{
\method{summary}{glasso.mix}(object,...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{object,...}{
an object with S3 class ``gm.select''. A list of the result from the function ``glassomix'' function.
  }
  
}

\details{
It gives a reduced summary of output from  ``glasso.mix''.
}

\value{
The details of the output components are as follows:
 \item{lambda }{The sequence of regularization parameters used.}
  \item{Pi}{ Mixture proportions  for each K across lambda.}
  \item{bestlambda.ebic}{Optimum lambda value based on EBIC for each K.}
  \item{besttheta.ebic}{The penalized precision matrix corresponding to the optimal EBIC for each value of K.}
  \item{n.par}{Total number of estimated parameters in the precision matrices corresponding to each 
             value of K at the various regularization parameters.}
}

\references{
Anani, Lotsi. and E. Wit (2012). High dimensional sparse Gaussian graphical model.
arXiv:1210.5371. \url{http://arxiv.org/abs/1210.5371}

Witten, Daniela M. and Friedman, Jerome H. and Simon, Noah (2011). New Insights and Faster 
Computations for the Graphical Lasso. \emph{Journal of Computational and Graphical Statistics}.
20(4), 892-900

Pan, Wei and Shen, Xiaotong (2007). Penalized Model-Based Clustering with Application to Variable Selection
\emph{J. Mach. Learn}. Res. 8(1532-4435), (1145--1164)

C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
  mclust Version 4 for R: Normal Mixture Modeling for Model-Based 
  Clustering, Classification, and Density Estimation. 
  Technical Report No. 597, Department of Statistics, University of Washington.
}
\author{
Anani Lotsi and Ernst Wit
}

\seealso{
  \code{\link{glasso.mix}}, 
  \code{\link{gm.select}}, 
  \code{\link{summary.gm.select}}, 
   \code{\link{gm.plot}} 
  }
\examples{
require(mvtnorm)
require(glasso)
require(igraph)
# Generate simulation data
p<-10   #number of graph nodes/variables             

Th1<-diag(p)
d=1
for (i in 1:p){
	for (j in 1:p){
		if ((abs(i-j)==d)|((i-j)/p ==d)|((j-i)/p==d)){
			Th1[i,j]=-0.4
		}
	}
}
	
Th2<-diag(p)
d=2
for (i in 1:p){
	for (j in 1:p){
		if ((abs(i-j)==d)|((i-j)/p ==d)|((j-i)/p==d)){
			Th2[i,j]=-0.4
		}
	}
}
	

n<-50   # number of replicates   
pi1 <- 0.15
z<-rbinom(n,1,pi1)

x<-NULL
for (i in 1:n){
	if (z[i]==1){
		x<-rbind(x,rmvnorm(1,rep(0,p),solve(Th1)))
	} else {
		x<-rbind(x,rmvnorm(1,rep(0,p),solve(Th2)))
	 }
}

ret=glasso.mix(x,K=NULL,lambda=NULL,em.iter=5,n.lambda=2,
                penalize.diagonal=TRUE,ebic.gamma=0.5,Kmax=3)
summary.glasso.mix(ret)
}

