\name{gm.plot}
\alias{gm.plot}
\title{
 Graphical plot  of the K Networks
}
\description{
Plots the K graphical precision matrices corresponding to the optimum EBIC.
}
\usage{
gm.plot(output)
}

\arguments{
  \item{output}{
It is a list which is the result of ``select.gm''  function.
}
  
}

\details{
It shows the graphical representation (dependencies ) of the p-variables in each cluster.
}


\references{
1.  Anani, Lotsi. and E. Wit (2013). High dimensional sparse Gaussian graphical model.
(http://arxiv.org/abs/1308.3381)


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
  \code{\link{select.gm}}, 
   \code{\link{summary.select.gm}}
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
graph=gm.plot(output)
}

