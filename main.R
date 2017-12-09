#E. Breza, A. Chandrasekhar, T. McCormick, and M. Pan
#Using aggregated relational data to feasibly identify network structure without network data.”


source('latent_surface_model.R')

# y - ARD data. n*K matrix, where n is the number of data points and K is the number of ARD questions.
# Can take NA as missing value 
# total.prop - the fraction of ties in the network that are made with members of group k, summed over K groups
# Can be approximated as the fraction of actors with feature k, summed over K groups
# muk.fix - fixed centers of subgroups on latent surface. m*p matrix, where m is the numebr of fixed centers, 
# and p is the number of dimension in latent space.
# n.iter, m.iter, and n.thin - parameters of the MCMC algorithm.
# is.sample - if we only have ARD data on a sample of the population
# distance.matrix - n.nonARD*n.ARD matrix. measure the distance in feature space between ARD and non-ARD nodes
# required if is.sample=TRUE
# Knn.K - parameter K in K nearest neighbors estimation
# output - a list of graphs: each simulated graph is n by n matrix, with ARD nodes before non-ARD nodes.

main=function(y,total.prop,muk.fix,n.iter=3000, m.iter=3, n.thin=10,is.sample=FALSE,distance.matrix=NULL,Knn.K=5,ls.dim=3){
  n=dim(y)[1]
  z.pos.init=generateRandomInitial(n,ls.dim)
  out=f.metro(y,total.prop=total.prop,n.iter=n.iter, m.iter=m.iter, n.thin=n.thin,z.pos.init=z.pos.init,muk.fix=muk.fix,ls.dim=ls.dim)
  posterior=getPosterior(out,n.iter,m.iter,n.thin,n)
  est.degrees=posterior$est.degrees
  est.eta=posterior$est.eta
  est.latent.pos=posterior$est.latent.pos
  est.gi=getGi(est.degrees,est.eta)
  if(is.sample){
    posteriorAll=getPosteriorAllnodes(distance.matrix,est.gi,est.latent.pos,Knn.K,ls.dim)
    est.gi.all=posteriorAll$est.gi.all
    est.latent.pos.all=posteriorAll$est.latent.pos.all
    g.sims=simulate.graph.all(est.degrees,est.eta,est.latent.pos,est.gi,est.gi.all,est.latent.pos.all,ls.dim)
  }else{
    g.sims=simulate.graph.all(est.degrees,est.eta,est.latent.pos,est.gi,est.gi,est.latent.pos,ls.dim)
  }
  return(g.sims)
}

generateRandomInitial=function(n,p){
  z=matrix(rnorm(n*p),nrow=n,ncol=p)
  z=sweep(z,MARGIN=1,1/sqrt(rowSums(z^2)),`*`)
  return (z)
}

getPosteriorAllnodes=function(distance.matrix,est.gi,est.latent.pos,Knn.K,ls.dim){
  n.ARD=dim(distance.matrix)[2]
  n.nonARD=dim(distance.matrix)[1]
  est.gi.all=NULL
  est.latent.pos.all=NULL
  for (ind in 1:dim(est.gi)[1]){
    g.ARD=est.gi[ind,]
    z.ARD=matrix(est.latent.pos[ind,],byrow=F,nrow=n.ARD,ncol=ls.dim)
    
    g.nonARD=NULL
    z.nonARD=NULL
    for (i in 1:n.nonARD){
      if(sort(distance.matrix[i,])[1]!=0){
        K.nn=order(distance.matrix[i,])[1:Knn.K]
        weights=(1/sort(distance.matrix[i,])[1:Knn.K])/sum(1/sort(distance.matrix[i,])[1:Knn.K])
        g.nonARD=c(g.nonARD,sum(g.ARD[K.nn]*weights))
        z.tmp=colSums(sweep(z.ARD[K.nn,],MARGIN=1,weights,'*'))
        z.nonARD=rbind(z.nonARD,z.tmp/sqrt(sum(z.tmp^2)))
      }else{
        g.nonARD=c(g.nonARD,g.ARD[order(distance.matrix[i,])[1]])
        z.nonARD=rbind(z.nonARD,z.ARD[order(distance.matrix[i,])[1]])
      }
    }
    g=c(g.ARD,g.nonARD)
    est.gi.all=rbind(est.gi.all,g)
    z=rbind(z.ARD,z.nonARD)
    est.latent.pos.all=rbind(est.latent.pos.all,c(z))
  }
  return(list(est.gi.all=est.gi.all,est.latent.pos.all=est.latent.pos.all))
}

getPosterior=function(out,n.iter,m.iter,n.thin,n){
  est.degrees=NULL
  est.eta=NULL
  est.latent.pos=NULL
  for(n.ind in (n.iter/n.thin/2+1):(n.iter/n.thin)){
    for (m.ind in 1:m.iter){
      est.degrees=rbind(est.degrees,out$sims[n.ind,m.ind,1:n])
      est.eta=c(est.eta,out$sims[n.ind,m.ind,][length(out$sims[n.ind,m.ind,])])
      est.latent.pos=rbind(est.latent.pos,c(out$sims.latent[n.ind,m.ind,1:n,]))
    }
  }
  return(list(est.degrees=est.degrees,est.eta=est.eta,est.latent.pos=est.latent.pos))
}

getGi=function(est.degrees,est.eta){
  gi.m=NULL
  for (ind in 1:length(est.eta)){
    nexp.gi=sqrt(sum(exp(est.degrees[ind,])) * cp.fcn(est.eta[ind])/cp.fcn(.000001))
    gi=exp(est.degrees[ind,])/nexp.gi* cp.fcn(est.eta[ind])/cp.fcn(.000001)
    gi.m=rbind(gi.m,gi)
  }
  return (log(gi.m))
}  

simulate.graph.all=function(est.degrees.ARD,est.eta,est.latent.pos.ARD,est.gi.ARD,est.gi,est.latent.pos,ls.dim){
  g.sims=list()
  n.ARD=dim(est.degrees.ARD)[2]
  n=dim(est.gi)[2]
  for (ind in 1:length(est.eta)){
    z=matrix(est.latent.pos[ind,],byrow=F,nrow=n,ncol=ls.dim)
    z.ARD=matrix(est.latent.pos.ARD[ind,],byrow=F,nrow=n.ARD,ncol=ls.dim)
    g.sims=c(g.sims,list(simulate.graph.once(z=z,g=est.gi[ind,],eta=est.eta[ind],d.ARD=est.degrees.ARD[ind,],z.ARD=z.ARD,g.ARD=est.gi.ARD[ind,])))
  }
  return(g.sims)
}


simulate.graph.once=function(z,g,eta,d.ARD,z.ARD,g.ARD){
  n.ARD=length(g.ARD)
  adjexp=matrix(NA,nrow=n.ARD,ncol=n.ARD)
  diag(adjexp)=0
  for(i in 1:(n.ARD-1)){
    for (j in (i+1):n.ARD){
      adjexp[i,j]=adjexp[j,i]=exp(g.ARD[i]+g.ARD[j]+eta*sum(z.ARD[i,]*z.ARD[j,]))
    }
  }
  const=sum(exp(d.ARD))/sum(adjexp)
  n=length(g)
  adj=matrix(NA,nrow=n,ncol=n)
  diag(adj)=0
  for(i in 1:(n-1)){
    for (j in (i+1):n){
      p.ij=exp(g[i]+g[j]+eta*sum(z[i,]*z[j,]))*const
      edge=rbinom(n=1,size=1,prob=min(p.ij,1))
      adj[i,j]=adj[j,i]=edge
    }
  }
  return(adj)
}



