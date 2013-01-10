\name{gm.plot}
\alias{gm.plot}
\title{
 Graphical plot (covariance structural ) of the K- Networks
}
\description{
This function is for ploting the K- graphical precision matrices corresponding to the optimum EBIC.
}
\usage{
gm.plot(output) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{output}{
a list which is the result of glasso.mix algorithm from the  ``glasso.mix.select'' function.
}

}
\details{
It shows the graphical representation (dependencies ) of the p-variables in each cluster.
}


\value{
The recovered K networks
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
  \code{\link{gm.select}},  
  \code{\link{summary.gm.select}}
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
graph=gm.plot(output)
}

