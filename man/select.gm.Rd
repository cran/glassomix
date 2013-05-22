\name{select.gm}
\alias{select.gm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{ High dimensional sparse Gaussian graphical mixture model selection

}
\description{
This function  selects the  optimal model according to Extended Bayesian Information Criterin (EBIC) for EM- algorithm for parameterized High dimensional 
sparse Gaussian graphical mixture models. The function etimates the optimun number of mixture components and the regularization parameter lambda.
}
\usage{
select.gm(ret)
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
  \item{n.cluster }{Optimal number of clusters or mixture components  K.}
  \item{eBIC }{All EBIC values. This is a matrix with row equals maximun number of mixtures minus one; (Kmax-1) and  column equals to n.lambda.}
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
1.  Anani, Lotsi. and E. Wit (2012). High dimensional sparse Gaussian graphical model.
(to be published)

2.  Witten, Daniela M. and Friedman, Jerome H. and Simon, Noah (2011). New Insights and Faster 
Computations for the Graphical Lasso. \emph{Journal of Computational and Graphical Statistics}.
20(4), 892-900

3.  Pan, Wei and Shen, Xiaotong (2007). Penalized Model-Based Clustering with Application to Variable Selection
\emph{J. Mach. Learn}. Res. 8(1532-4435).

4.  C. Fraley, A. E. Raftery, T. B. Murphy and L. Scrucca (2012).
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
   \code{\link{summary.select.gm}},
  \code{\link{gm.plot}}
 
  }
    
\examples{
#generate data 
p<-10          ## Number of observations
n<-100         ## Number of nodes

L1 = huge.generator(n = n, d = p, vis = FALSE, graph = "random")
L2 = huge.generator(n = n, d = p, vis = FALSE, graph = "random")
Th1<-L1$sigma  ## Precision matrix from graph1
Th2<-L2$sigma   ## Precision matrix from graph2
pi1 <- 0.5
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
               
output=select.gm(ret)
}

