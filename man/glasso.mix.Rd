\name{glasso.mix}
\alias{glasso.mix}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 sparse Gaussian undirected graphical mixture model estimation.
}
\description{
The main function perfoming the inference via EM algorithm. This function for each value 
of K, estimate the responsibility matrices ( \eqn{n \times K}{n x K}) at the E-step 
and then given these probabilities, estimate the precision matrices at the M-step 
via glasso.
}
\usage{
glasso.mix(data,K=NULL,lambda=NULL,em.iter,n.lambda,
           penalize.diagonal=TRUE,ebic.gamma=0.5,Kmax)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
 ( \eqn{n \times p}{n x p}),  
        rows = n, number of observation, 
        columns = p, number of graph nodes/variables).
}
  \item{K}{
 A sequence of integers denoting the numbers of mixture components (clusters) for which  EBIC is to be calculated. 
 }
  \item{lambda}{
(Non-negative) regularization parameter for glasso. lambda=0 means no regularization. It could be a scalar  or a vector. 
}
  \item{em.iter}{
The maximun number of EM iteration. 
}
  \item{n.lambda}{
The length of the tuning parameter lambda.
}
  \item{penalize.diagonal}{
Should diagonal of precision matrice be penalized? Dafault is FALSE.
}
  \item{ebic.gamma}{
The Extended Bayesian Information Criteria paremeter, usually  ebic.gamma is between 0 and 1. 
}
  \item{Kmax}{
The maximum number of K. 
}
}
\details{
 Implements the EM algorithm for a parameterized Gaussian graphical mixture models accross K for each of the regularization parameters.
 
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
1.  Anani, Lotsi. and E. Wit (2013). High dimensional sparse Gaussian graphical model.
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
  \code{\link{summary.glasso.mix}},
  \code{\link{select.gm}}, 
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

## Output from the main function
ret=glasso.mix(x,K=NULL,lambda=NULL,em.iter=5,n.lambda=2,
               penalize.diagonal=TRUE,ebic.gamma=0.5,Kmax=3)
}




