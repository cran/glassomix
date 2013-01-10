\name{glasso.mix}
\alias{glasso.mix}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
EM algorithm  for a parameterized sparse Gaussian graphical mixture model.
}
\description{
This function perfoms the inference via EM algorithm. This function for each value 
of K, estimate the responsibility matrices ( \eqn{n \times K}{n x K}) at the E-step 
and then given these probabilities,  res-estimate the precision matrices at the M-step 
via glasso.
}
\usage{
glasso.mix(data,K=NULL,lambda=NULL,em.iter,n.lambda,
                penalize.diagonal=TRUE,ebic.gamma=0.5,Kmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
It is a matrix of the data ( \eqn{n \times p}{n x p} matrix) 
        rows = n, number of replicates, 
        columns = p, number of graph nodes/variables).
}
  \item{K}{
Number of mixture components. A sequence of integers denoting the numbers of mixture components (clusters) for which the EBIC is to be calculated. Generaly it is a sequence of 2, 3,...whose lenght should equal the length of lambda+1
}
  \item{lambda}{
(Non-negative) regularization parameter for glasso. lambda=0 means no regularization. Can be a scalar  or a vector. 
}
  \item{em.iter}{
The maximun number of EM iteration. Usually the convergence is faster.
}
  \item{n.lambda}{
The length of the tuning parameter lambda.
}
  \item{penalize.diagonal}{
Should diagonal of precision matrice be penalized? Dafault FALSE.
}
  \item{ebic.gamma}{
The Extended Bayesian Information Criteria paremeter, usually  ebic.gamma is between 0 and 1. 
}
  \item{Kmax}{
The maximum number of K 
}
}
\details{
 Implements the EM algorithm for a parameterized Gaussian graphical mixture model accross K for each of the regularization parameters,
  starting with the expectation step.
}

\value{
The details of the output components are as follows:
  \item{res }{
    A list with the following components:
      \describe{ 
         \item{\code{loglik}}{
             A vector value of un-penalized log-likelihood for each value of K.
         }
               
         \item{\code{naiveloglik}}{
             A vector value of  naive log-likelihood extracted from glasso for each value of K.
         }
          
        \item{\code{n.par}}{
             Total number of estimated parameters in each the precision matrices corresponding to each 
             value of K at the various regularization parameters.
         }
                
        \item{\code{bestlambda.ebic}}{
             Optimal tuning parameter corresponding to each 
             value of K.
         }
         
        \item{\code{besttheta.ebic}}{
             The penalized precision matrix corresponding to the optimal EBIC for each value of K.
         }
         
         \item{\code{bestpi.ebic}}{
             The mixture proportion corresponding to the optimal EBIC for each value of K.
         }  
        \item{\code{Theta_Pen}}{
             All the individual penalized precision matrices corresponding to  each value of K 
             at the various regularization parameters.
         }
             
        \item{\code{Theta_NonPen}}{
             All the individual non-penalized precision matrices corresponding to  each value of K 
             at the various regularization parameters.
         }
            
        \item{\code{pi.ind}}{
             Responsibility matrices ( \eqn{n \times K}{n x K}) corresponding to  each value of K for the various regularization parameters.
             It can also be seen as vector of probabilities (\eqn{w^{(k)}_{i1},...,w^{(K)}_{iK}}{w^{(k)}_{i1},...,w^{(K)}_{iK}}) 
             of individual \eqn{i}{i} belonging to the \eqn{k}{k} classes at penalty  \eqn{\lambda}{lambda}.
         }
            
        \item{\code{pi}}{
             K- Mixing coefficients  for the various regularization parameters.
         }
             
        \item{\code{EBIC}}{
             All EBIC values for each value of K.
         }   
         
                   
    }          
  }
 \item{lambda }{The sequence of regularization parameters used.}
 \item{Kmax }{The maximun number of mixture components.}
  \item{n.lambda }{The length of the tuning parameter lambda.}
   \item{data}{The data matrix.}

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
 \code{\link{summary.glasso.mix}},
  \code{\link{gm.select}}, 
   \code{\link{summary.gm.select}},
  \code{\link{gm.plot}}
   
  }

\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
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

}




