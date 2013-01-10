\name{glassomix-package}
\alias{glassomix-package}
\alias{glassomix}
\docType{package}
\title{
High dimensional sparse Gaussian graphical mixture model
}
\description{
The R package ``glassomix'' is a statistical tool for performing model-based soft 
clustering. It provides functions for parameter estimation via the EM algorithm 
for Gaussian graphical mixture models  in high dimensional setting. 
The main function is ``glasso.mix'' upon which a model selection is performed. 
The package estimates the optimum number of mixture components (K) and the tuning 
parameter (lambda) based on the Extended Bayesian Information Criteria (EBIC) via 
``gm.select'' function. The  graphical structural of the K-networks are also 
plotted through the function ``gm.plot''. 
}

\details{
This package provides an implementation of the procedures described in Lotsi and Wit (2013). 
The main function is ``glasso.mix'' which is a high dimensional algorithm for 
Gaussian graphical mixture model based on Extended  Bayesian Information Criterion.
The graphical structural of the K-subgroup of population of individuals is estimated 
and plotted. The package basically seeks to recover the optmum number of mixture components (K)
and then discern and discover the dependicies within the subgroup of individuals  
with no prior information about which cluster an individual propely belongs.

Functions:
  \preformatted{
glasso.mix            Inference via EM algorithm. This function for each value of K, 
                      estimates the responsibility matrices (n by K) at the E-step and
                      then given these probabilities, res-estimate the precision matrices 
                      at the M-step via glasso.  
                      
summary.glasso.mix    Summary of the result according to  fuction ``glasso.mix  ''  
                    
gm.select             This function performs the model selection and selects the 
                      best graphical models based on EBIC.
                      
summary.gm.select     Summary of the result according to ``gm.selec''
 
gm.plot               Plot the resultant graphical structural (covariance structural ) of the
                      K- subgroup of population of individuals. 


  }
}
\author{
Anani Lotsi <a.lotsi@rug.nl> and Ernst Wit
}
\references{
Anani, Lotsi. and E. Wit (2012). High dimensional sparse Gaussian graphical model.
arXiv:1210.5371. \url{http://arxiv.org/abs/1210.5371}

Witten, Daniela M. and Friedman, Jerome H. and Simon, Noah (2011). New Insights and Faster 
Computations for the Graphical Lasso. \emph{Journal of Computational and Graphical Statistics}.
20(4), 892-900

Pan, Wei and Shen, Xiaotong (2007). Penalized Model-Based Clustering with Application to Variable Selection
\emph{J. Mach. Learn. Res}. 8(1532-4435), (1145--1164)
}
\keyword{package }
