\name{gm.select}
\alias{gm.select}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ High dimensional sparse Gaussian graphical mixture model

}
\description{
This function performs the model selection and select the  optimal model according to Extended Bayesian Information Criterin (EBIC) for EM- algorithm for parameterized High dimensional 
sparse Gaussian graphical mixture models. The function etimates the optimun number of mixture components i.e number of clusters and the regularization parameter lambda.
}
\usage{
gm.select(ret)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ret}{
It is a list which is the result of glasso.mix algorithm from  'glasso.mix' function.
}
}
\details{
Implements the model selection clustering through a model selection based on  the EBIC for a parameterized 
Gaussian graphical mixture model accross K for each of the regularization parameters.
}
\value{
The details of the output components are as follows:
  \item{n.cluster }{Optimal number of clusters mixture components  K.}
  \item{eBIC }{All EBIC values. This is a matrix with row equals maximun number of mixtures minus one; (Kmax-1) and  column equals n.lambda.}
  \item{lambda_eBIC}{Optimum lambda value based on EBIC.}
  \item{Th.Pen}{Optimum penalized K precision matrices.}
  \item{Th.NPen}{Optimum non-penalized K precision matrices.}
  \item{Pi_ind}{Optimum  responsibility matrices ( \eqn{n \times K}{n x K}) corresponding to the soft-K-means clustering.}
  \item{Pi}{Optimum (K) mixture proportions  based on EBIC criterion.}
  \item{clusters}{( \eqn{n \times 1}{n x 1}) vector containing the indices of the clusters where the data points are assigned to.}
  \item{Pen_LogLik}{The un-penalized loglikelihood corresponding to the optimal EBIC.}
  \item{NPen_LogLik}{The naive un-penalized loglikelihood corresponding to the optimal EBIC.}
  \item{lambda}{The sequence of regularization parameters used.}
   
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
    \code{\link{summary.glasso.mix}},
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
output=gm.select(ret)
}

