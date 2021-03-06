\name{glassomix-package}
\alias{glassomix-package}
\alias{glassomix}
\docType{package}
\title{
High dimensional undirected Graphical Mixture Models selection
}
\description{
A package for high  dimensional Undirected Graphical Mixture Models selection. 
}

\details{
This package provides an implementation of the procedures described in Lotsi and Wit (2013)
(http://arxiv.org/abs/1308.3381). The main function is ``glasso.mix''. 
This function performs the graph estimation using glasso and a model selction 
is performed based on Extended  Bayesian Information Criterion  through the function
``select.gm''. The graphical structural of the K-subgroup of population  of individuals
is estimated and plotted via the function ``gm.plot''. 

Functions:
  \preformatted{
glasso.mix            Inference via EM algorithm based on glasso. 
                      
summary.glasso.mix    Summary of the result according to  fuction ``glasso.mix ''  
                    
select.gm             This function performs the model selection and selects the 
                      best graphical models based on EBIC.
                      
summary.select.gm     Summary of the result according to ``select.gm ''
 
gm.plot               Plot the resultant graphical structural  of the
                      K- subgroup of population of individuals. 
}
}
\author{
Anani Lotsi <a.lotsi@rug.nl> and Ernst Wit
}
\references{
Anani, Lotsi. and E. Wit (2013). Sparse  Gaussian Graphical Mixture Model.
(http://arxiv.org/abs/1308.3381)
}
\keyword{package }
