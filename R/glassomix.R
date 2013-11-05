
## All functions for the "glasso.mix" package
####################################################
# Inference via EM algorithm
####################################################
glasso.mix <- function(data,K=NULL,lambda=NULL,em.iter,n.lambda,penalize.diagonal=TRUE,ebic.gamma=0.5,Kmax){
	# data = data matrix (rows = n, number of replicates
	#                  cols = p, number of graph nodes/variables)
	# K = number of mixtures
	# em.iter = number of EM iterations
	#
   start.time <- Sys.time()
   if (is.matrix(data) == F){
    stop("Data should be a matrix or dataframe")
   }
   if (any(is.na(data))) stop("Data should contain no missing data")
   if (is.null(K)){
    K<-seq(1,Kmax)
  }
  res<-list()
  for(K in 1:Kmax){
    n<-dim(data)[1]
	  p<-dim(data)[2]
    # Initial values
	  pi.hat<-rep(1/K,K)
	  Th.hat<-array(NA,c(p,p,K))
    ThNP.hat<-array(NA,c(p,p,K))
	  Sigma.hat<-array(NA,c(p,p,K))
	   for (k in 1:K){
         if(n<p){
              Sigma.hat[,,k]<-diag(jitter(apply(data,2,var)))
            } else{
              S<-cov(data)
              if(n>p)Sigma.hat[,,k]<-rWishart(1,2*p,S)[,,1]
            }

		Th.hat[,,k]<-solve(Sigma.hat[,,k])
	}
	pi.ind<-matrix(NA,ncol=K,nrow=n)
  penloglik<-matrix(NA,nrow=em.iter,ncol=K)
	naiveNPLik<-matrix(NA,nrow=em.iter,ncol=K)
  if (is.null(lambda)){
    s2<-sort(apply(data,2,var),T)
    lambda<-seq(0.01,sqrt(s2[1]*s2[2])/2,length=n.lambda)
  }
  Rho<-ebic.gamma
	Th.lambda<-list(array(NA,c(p,p,K)))
	ThNP.lambda<-list(array(NA,c(p,p,K)))
	pi=NULL
	pi.i=NULL
	loglik.lambda<-vector("list", length(lambda))
  NPLik.lambda<-vector("list", length(lambda))
  naiveNPLik.lambda<-vector("list", length(lambda))
  param<-vector()
  EBIC = vector()
  naiveBIC = vector()
  naiveEBIC = vector()

   ##  for search of initial values

		      # E step
  for (k in 1:K){
	   pi.ind[,k]<-dmvnorm(data,rep(0,p),Sigma.hat[,,k])*pi.hat[k]
  }
  pi.ind<-pi.ind/apply(pi.ind,1,sum)
  pi.hat<-apply(pi.ind,2,mean)
  w<-t(t(pi.ind)/apply(pi.ind,2,sum))
        # M step
  for (k in 1:K){
    Sw<-t(data)%*%diag(w[,k])%*%data
    Thk.gl <-glasso(Sw,rho=10,penalize.diagonal=penalize.diagonal)     		# M step
    Th.hat[,,k] <- Thk.gl$wi
    Sigma.hat[,,k] <- Thk.gl$w
  }
                 # End of search  for initial values

	for(l in 1:length(lambda)){      ## Iterations accross lambda

     ## Begin the EM  (Penalized EM)
      for (i in 1:em.iter){
		      # E step
		      for (k in 1:K){
		      	 pi.ind[,k]<-dmvnorm(data,rep(0,p),Sigma.hat[,,k])*pi.hat[k]
		      }
          pi.ind<-pi.ind/apply(pi.ind,1,sum)
		      pi.hat<-apply(pi.ind,2,mean)

          w<-t(t(pi.ind)/apply(pi.ind,2,sum))

          # M step
		      for (k in 1:K){
			       Sw<-t(data)%*%diag(w[,k])%*%data
			       Thk.gl <-glasso(Sw,(2*lambda[l]/max(sum(pi.ind[,k]),0.001)),penalize.diagonal=penalize.diagonal)     		# M step
			       penloglik[i,k]<-Thk.gl$loglik
			       Th.hat[,,k] <- Thk.gl$wi
			       Sigma.hat[,,k] <- Thk.gl$w
             naiveNPLik[i,k] = (1/2)*(determinant(Thk.gl$wi)$modulus[1] - sum(Thk.gl$wi*t(Sw))) ##Compute un-penalized likelihood
          }
        }
      ## Begin the EM (Unpenalized EM)
      for (i in 1:em.iter){
		      # E step
		      for (k in 1:K){
		     	  pi.ind[,k]<-dmvnorm(data,rep(0,p),Sigma.hat[,,k])*pi.hat[k]
		      }
          pi.ind<-pi.ind/apply(pi.ind,1,sum)
		      pi.hat<-apply(pi.ind,2,mean)

          w<-t(t(pi.ind)/apply(pi.ind,2,sum))
          NPLik<-matrix(nrow=K,ncol=n)
		      for (k in 1:K){
	           Sw<-t(data)%*%diag(w[,k])%*%data
              if (sum(Th.hat[,,k]==0)==0){
			           Thk.gl <-glasso(Sw,0.001,penalize.diagonal=penalize.diagonal)
              }else{
			           Thk.gl <-glasso(Sw,0.001,penalize.diagonal=penalize.diagonal,zero=which(Th.hat[,,k]==0,arr.ind=T))
			       }
			       Sigma.hat[,,k] <- Thk.gl$w
             ThNP.hat[,,k] <- Thk.gl$wi   		# M step
			       wx<-data*sqrt(w[,k])
             lx<-apply(wx,1,function(xx,Th){matrix(xx,nrow=1)%*%Th%*%matrix(xx,ncol=1)},Th=Thk.gl$wi)
             ldet<- determinant(Thk.gl$wi)$modulus[1]
	           NPLik[k,]<- (1/2)*(ldet - lx)
          }
       }
        rownames(Th.hat) <- colnames(data)
        colnames(Th.hat) <- colnames(data)
		    Th.lambda[[l]]<-Th.hat
		    ThNP.lambda[[l]]<-ThNP.hat
		    pi[[l]]=pi.hat
		    pi.i[[l]]=pi.ind
		    NPLik.lambda[[l]]<-sum(log(apply(pi.hat*exp(NPLik),2,sum)))
        naiveNPLik.lambda[[l]]<-apply(naiveNPLik,1,sum)[em.iter]
        param[l]<-  sum(Th.lambda[[l]]!=0)
        EBIC[l] = -2*NPLik.lambda[[l]] + log(n)*(param[l]/2)+4*param[l]*Rho*log(p)     ## Extended BIC
        naiveEBIC[l] = -2*naiveNPLik.lambda[[l]] +log(n)*(param[l]/2)+4*param[l]*Rho*log(p)     ## Extended BIC


      }    # End of lambda
      besttheta.ebic<- Th.lambda[[which.min(EBIC)]]
      bestlambda.ebic<-lambda[which.min(EBIC)]
      bestpi.ebic<- pi[[which.min(EBIC)]]
      res[[K]]<-list(loglik=NPLik.lambda,naiveloglik=naiveNPLik.lambda,n.par=param, bestlambda.ebic= bestlambda.ebic,
                besttheta.ebic=besttheta.ebic,bestpi.ebic=bestpi.ebic,Theta.Pen=Th.lambda,Theta.NonPen=ThNP.lambda,pi.ind=pi.i,pi=pi,EBIC=EBIC)
      ret<-list(res=res,lambda=lambda,Kmax=Kmax,n.lambda=n.lambda,data=data)
   } ############# END of K
   class(ret) <- "glasso.mix"
return(ret)
}


######################### Summary function based on  glasso.mix. Reduced summary from glasso.mix   #####################################

summary.glasso.mix<-function(object,...){
  pi<-NULL
  bestlambda.ebic<-NULL
  besttheta.ebic<-NULL
  n.par<-NULL
  Kmax<-object$Kmax
  lambda<-object$lambda
      for(k in 1:Kmax){
          pi[[k]]<-object$res[[k]]$pi
          bestlambda.ebic[[k]]<-object$res[[k]]$bestlambda.ebic
          besttheta.ebic[[k]]<-object$res[[k]]$besttheta.ebic
          n.par[[k]]<-object$res[[k]]$n.par
      }
  return.list<-list(lambda=lambda,pi=pi,bestlambda.ebic=bestlambda.ebic,besttheta.ebic=besttheta.ebic,n.par=n.par)
  return(return.list)

}
###############################################################################################################

###############################################################################################################

  #The function below does the model selection (hat{K},hat{lambda})=argmin_{K,lambda}(EBIC) and output the followings:
  #(hat{K},hat{lambda}):  The pairs that minimizes the EBIC
  #The "K" penalized precision matrices:  "Pen_LogLik"
  #The "K" Un-penalized precision matrices:  "NPen_LogLik"
  #bestpi.ind: An (n by K) soft clustes
  #clusters: The individial clustering

###################################################################################################################
select.gm<-function(ret){
    cl=NULL   ## the class of the elelements of the sample
    #ret=glasso.mix(data,K=NULL,lambda=NULL,em.iter,n.lambda,penalize.diagonal=TRUE,ebic.gamma=0.5,Kmax)
    lambda<-ret$lambda
    Kmax<-ret$Kmax
    n.lambda<-ret$n.lambda
    data<-ret$data
    EBICklambda = matrix(0, nrow=Kmax-1,ncol=n.lambda)
    colnames(EBICklambda)=round(lambda,2)
    rownames(EBICklambda)=seq(2,Kmax)
    for(K in seq(2,Kmax)){EBICklambda[K-1,]<-ret$res[[K]]$EBIC }
    bestklambda<-which(EBICklambda == min(EBICklambda), arr.ind = TRUE)
    k_index<-bestklambda[1,][1]
    n.cluster=seq(2,Kmax)[k_index]
    lambda_index<-bestklambda[1,][2]
    lambda_eBIC<-lambda[lambda_index]
    Pi_ind= ret$res[[n.cluster]]$pi.ind[[lambda_index]]
    Pi<-apply(Pi_ind,2,mean)
    for(a in 1:dim(data)[1]){
        cl[a]<-which.max(Pi_ind[a,])
   }

  output<-list(n.cluster=n.cluster,eBIC=EBICklambda, lambda.eBIC=lambda_eBIC,Th.Pen= ret$res[[n.cluster]]$Theta.Pen[[lambda_index]],Th.NPen= ret$res[[n.cluster]]$Theta.NonPen[[lambda_index]],
  Pi.ind=ret$res[[n.cluster]]$pi.ind[[lambda_index]],Pi=Pi,clusters=cl,Pen.LogLik=ret$res[[n.cluster]]$loglik[lambda_index], NPen.LogLik=ret$res[[n.cluster]]$naiveloglik[lambda_index],lambda=lambda)

  class(output) <- "select.gm"
  return(output)

}


################     summary of the result according to gm.select ###############
summary.select.gm<-function(object,...){
  mix_comp<-object$n.cluster
  lambda_eBIC<-object$lambda_eBIC
  clustering<-object$clusters
  mix_prop<-object$Pi
  return.list<-list(mix.comp= mix_comp,lambda.eBIC=lambda_eBIC, clustering= clustering,mix.prop=mix_prop)
  return(return.list)
}

#######################################################################

######## Graphical plot ####################
#This function plots the "K" networks
##########################################
gm.plot<-function(output) {
    net=output$Th.Pen
    K=output$n.cluster
    par(mfrow = c(1, K))
    for(k in 1:K){
      g1<-1*(abs(net[,,k])>0.01)
      g1 <- graph.adjacency(g1,mode = "undirected" )
      plot.igraph(g1, layout = layout.circle, main = "network")

    }
}

##########   END ###############################################################################



