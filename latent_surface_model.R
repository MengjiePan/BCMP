#E. Breza, A. Chandrasekhar, T. McCormick, and M. Pan
#Using aggregated relational data to feasibly identify network structure without network data.”

#Priors used in this version of code: uniform flat priors for a.mu, a.s, b.mu, and b.s; 
#gamma(0.5,0.5) for v.eta; gamma(5,0.1) for v.etak. Experiments on how different priors affect the accuracy
#of estimation can be found in appendix D of the paper.

rot.mat<-function(r.v){
  phi=r.v[1]
  theta=r.v[2]
  psi=r.v[3]
  rot.tmp<-matrix(NA,3,3)
  rot.tmp[1,]<-c(cos(theta)*cos(psi),((-cos(phi)*sin(psi))+(sin(phi)*sin(theta)*cos(psi))),((sin(phi)*sin(theta))+(cos(phi)*sin(theta)*cos(psi))))
  rot.tmp[2,]<-c(cos(theta)*sin(psi),((cos(phi)*cos(psi))+(sin(phi)*sin(theta)*sin(psi))),((-sin(phi)*cos(psi))+(cos(phi)*sin(theta)*sin(psi))))
  rot.tmp[3,]<-c(-sin(theta),sin(phi)*cos(theta),cos(phi)*cos(theta))
  return(rot.tmp)
}

sphere.update<-function(u.mean,u.k,m=3,rot=F,rot.vec=NULL){
  #step0
  b=((-2*u.k)+((4*(u.k)^2)+(m-1)^2)^(1/2))/(m-1)
  x0=(1-b)/(1+b)
  c=(u.k*x0)+((m-1)*log(1-(x0)^2))
  #step1
  z.tmp.u<-rbeta(length(u.k),(m-1)/2,(m-1)/2)
  u.tmp.u<-runif(length(u.k),0,1)
  w.u<-(1-((1+b)*z.tmp.u))/(1-((1-b)*z.tmp.u))
  #step2
  iter=0
  rejump.ind<-which(log(u.tmp.u)>((u.k*w.u)+((m-1)*log(1-(x0*w.u)))))
  while(sum(log(u.tmp.u)>((u.k*w.u)+((m-1)*log(1-(x0*w.u)))))>0 & iter<200){
    z.tmp.u.tmp<-rbeta(length(rejump.ind),(m-1)/2,(m-1)/2)
    u.tmp.u.tmp<-runif(length(rejump.ind),0,1)
    w.u.tmp<-(1-((1+b[rejump.ind])*z.tmp.u.tmp))/(1-((1-b[rejump.ind])*z.tmp.u.tmp))
    iter=iter+1
    if(iter>100){print(iter)}
    w.u[rejump.ind]<-w.u.tmp
    u.tmp.u[rejump.ind]<-u.tmp.u.tmp
    z.tmp.u[rejump.ind]<-z.tmp.u.tmp
    rejump.ind<-which(log(u.tmp.u)>((u.k*w.u)+((m-1)*log(1-(x0*w.u)))))
  }
  #step3
  #unit vector
  v.tmp.u=matrix(runif((length(u.k)*(m-1)),-100,100),length(u.k),m-1)
  v.tmp.u=v.tmp.u/sqrt(rowSums(v.tmp.u^2))
  xt=cbind(((1-w.u^2)^(1/2))*v.tmp.u,w.u)
  if(rot==T){
    if(is.null(rot.vec)==F){
      return(rot.mat(rot.vec)%*%t(xt))}
    if(is.null(rot.vec)==T){
      den=c(u.mean[1],u.mean[2])
      den=sqrt(t(den)%*%(den))
      r.tmp=c(0,acos(u.mean[3]),acos(u.mean[1]/den))
      if(u.mean[2]<0){r.tmp=c(0,acos(u.mean[3]),-acos(u.mean[1]/den))}
      return(rot.mat(r.tmp)%*%t(xt))
    }}
  if(rot==F){return(xt)}
}
cp.fcn<-function(kap,p=3){
  top<-kap^((p/2)-1)
  bot<-((2*pi)^(p/2))*besselI(kap,((p/2)-1))
  out<-top/bot
  return(out)
}

f.datalik<-function(data, v.alpha, v.beta, v.eta, v.etak, v.z, v.muk){
  data.n<-nrow(data)
  data.j<-ncol(data)
  
  v.theta<-v.z%*%t(v.muk)
  v.null<-exp(outer(v.alpha, v.beta, "+"))
  v.pop<-cp.fcn(v.eta)
  v.distr<-cp.fcn(v.etak)/(cp.fcn(.000001)*cp.fcn(sqrt((v.eta^2)+(v.etak^2)+(2*v.eta*v.etak*cos(v.theta)))))
  v.lambda<-v.null*v.pop*v.distr

  prob.y<-(-v.lambda)+(data*log(v.lambda))-log(factorial(data))
  return(prob.y)
}


f.update<-function(data, v.alpha, v.beta, v.eta, v.etak, v.z, v.muk, a.mu, a.s, b.mu, b.s, ad.jump, bd.jump,eta.jump,etak.jump,muk.jump,z.jump,k.mu,k.s,vjind, total.prop,muk.fix,ls.dim){
  data.n<-nrow(data)
  data.j<-ncol(data)
  
  ### updating alpha
  a.star<-v.alpha+rnorm(length(v.alpha),0, sqrt(ad.jump))
  a.old<-v.alpha
  lik.old<-f.datalik(data, a.old, v.beta, v.eta, v.etak, v.z, v.muk)
  lik.star<-f.datalik(data, a.star, v.beta,v.eta, v.etak, v.z, v.muk)
  prob.alpha.diff<-rowSums(lik.star, na.rm=T)+log(dnorm(a.star, a.mu, sqrt(a.s)))-rowSums(lik.old, na.rm=T)-log(dnorm(a.old, a.mu, sqrt(a.s)))
  alpha.valid<-(!is.infinite(prob.alpha.diff))&(!(is.na(prob.alpha.diff)))
  jump.new<-rep(0, length(v.alpha))   
  jump.new[alpha.valid]<-rbinom(sum(alpha.valid), 1, exp(pmin(prob.alpha.diff[alpha.valid], 0)))
  
  v.alpha<-ifelse(jump.new, a.star, a.old)
  alpha.p<-rep(0, length(v.alpha))
  alpha.p[alpha.valid]<- exp(pmin(prob.alpha.diff[alpha.valid], 0))
  
  ### updating a.mu
  a.mu<-rnorm(1, mean(v.alpha, na.rm=T), sqrt(a.s/data.n))
  ### updating a.s
  a.s=sum((v.alpha-a.mu)^2, na.rm=T)/rchisq(1, data.n-1)
  ### updating beta
  b.new<-v.beta+rnorm(length(v.beta),0, sqrt(bd.jump))
  lik.old<-f.datalik(data, v.alpha, v.beta, v.eta, v.etak, v.z, v.muk)
  b.old<-v.beta
  v.beta<-b.new
  lik.new<-f.datalik(data, v.alpha, v.beta, v.eta, v.etak, v.z, v.muk)
  prob.beta.diff<-colSums(lik.new, na.rm=T)+log(dnorm(b.new, b.mu, sqrt(b.s)))-colSums(lik.old, na.rm=T)-log(dnorm(b.old, b.mu, sqrt(b.s)))
  beta.valid<-(!is.infinite(prob.beta.diff))&(!is.na(prob.beta.diff))
  jump.new<-rep(0, length(v.beta))
  beta.p<-rep(0, length(v.beta))
  v.beta<-b.old
  if(sum(beta.valid)>0){
    jump.new[beta.valid]<-rbinom(sum(beta.valid), 1, exp(pmin(prob.beta.diff[beta.valid], 0)))
    v.beta<-jump.new*b.new+(1-jump.new)*b.old
    beta.p[beta.valid]<- exp(pmin(prob.beta.diff[beta.valid], 0))}
  
  ### updating b.mu
  b.mu<-rnorm(1, mean(v.beta), sqrt(b.s/data.j))
  ### updating b.s
  b.s<-sum((v.beta-b.mu)^2)/rchisq(1, data.j-1)
  ### renormalize
  const<-log(sum(exp(v.beta[1:8]))/total.prop)
  v.alpha<-v.alpha+const
  a.mu<-a.mu+const
  v.beta<-v.beta-const
  b.mu<-b.mu-const
  
  ### updating eta
  eta.star<-abs(v.eta+rnorm(1,0, sqrt(eta.jump)))
  eta.old<-v.eta
  lik.old<-f.datalik(data, v.alpha, v.beta, eta.old, v.etak, v.z, v.muk)
  lik.star<-f.datalik(data, v.alpha, v.beta,eta.star, v.etak, v.z, v.muk)
  prob.eta.diff<-sum(lik.star, na.rm=T)+log(dgamma(eta.star, .5,.5))-sum(lik.old, na.rm=T)-log(dgamma(eta.old, .5,.5))
  eta.valid<-(!is.infinite(prob.eta.diff))&(!(is.na(prob.eta.diff)))
  jump.new<-0
  eta.p<-0
  v.eta<-eta.old
  if(eta.valid==T){   
    jump.new<-rbinom(1, 1, exp(pmin(prob.eta.diff, 0)))
    v.eta<-ifelse(jump.new, eta.star, eta.old)
    eta.p<- exp(min(prob.eta.diff[eta.valid], 0))}
  
  ## updating etak
  etak.new<-abs(v.etak+rnorm(length(v.etak),0, sqrt(20)))
  etak.new[etak.new>250]<-abs(250-abs(etak.new[etak.new>250]-250))
  lik.old<-f.datalik(data, v.alpha, v.beta, v.eta, v.etak, v.z, v.muk)
  etak.old<-v.etak
  v.etak<-etak.new
  lik.new<-f.datalik(data, v.alpha, v.beta, v.eta, v.etak, v.z, v.muk)
  prob.etak.diff<-colSums(lik.new, na.rm=T)+log(dgamma(etak.new,5,.1))-colSums(lik.old, na.rm=T)-log(dgamma(etak.old,5,.1))
  etak.valid<-(!is.infinite(prob.etak.diff))&(!is.na(prob.etak.diff))&(etak.new<1450)
  jump.new<-rep(0, length(v.etak)) 
  etak.p<-rep(0, length(v.etak))
  v.etak<-etak.old
  if(sum(etak.valid)>0){   
    jump.new[etak.valid]<-rbinom(sum(etak.valid), 1, exp(pmin(prob.etak.diff[etak.valid], 0)))
    v.etak<-jump.new*etak.new+(1-jump.new)*etak.old
    etak.p[etak.valid]<- exp(pmin(prob.etak.diff[etak.valid], 0))}
  
  k.mu<-rnorm(1, mean(log(v.etak)), sqrt(k.s/data.j))
  k.s<-sum((log(v.etak)-k.mu)^2)/rchisq(1, data.j-1)

  ### updating z
  z.star<-array(dim=dim(v.z))
  for(mm in 1:nrow(z.star)){
    z.star[mm,]<-sphere.update(u.mean=v.z[mm,],u.k=z.jump[mm],m=ls.dim,rot=T)}
  z.old<-v.z
  lik.old<-f.datalik(data, v.alpha, v.beta, v.eta, v.etak, v.z, v.muk)
  lik.star<-f.datalik(data, v.alpha, v.beta,v.eta, v.etak, z.star, v.muk)
  prob.z.diff<-rowSums(lik.star, na.rm=T)-rowSums(lik.old, na.rm=T)
  z.valid<-(!is.infinite(prob.z.diff))&(!(is.na(prob.z.diff)))
  jump.new<-rep(0, nrow(v.z))   
  jump.new[z.valid]<-rbinom(sum(z.valid), 1, exp(pmin(prob.z.diff[z.valid], 0)))
  v.z[which(jump.new==1),]<-z.star[which(jump.new==1),]
  z.p<-rep(0, nrow(v.z))
  z.p[z.valid]<- exp(pmin(prob.z.diff[z.valid], 0))

  ### updating mukk
  mumean=colMeans(muk.fix)/sqrt(colMeans(muk.fix)%*%colMeans(muk.fix))
  muk.old<-v.muk
  muk.new<-array(dim=dim(v.muk))
  for(m in 1:nrow(muk.new)){
    muk.new[m,]<-rmovMF(n=1,theta=((.1*mumean)+(muk.old[m,]*v.etak[m])),alpha=1)}
  muk.p<-rep(0, nrow(v.muk))
  v.muk=muk.new
  v.muk[muk.fix.ind[1],]<-muk.fix[1,]
  v.muk[muk.fix.ind[2],]<-muk.fix[2,]
	v.muk[muk.fix.ind[3],]<-muk.fix[3,]
	v.muk[muk.fix.ind[4],]<-muk.fix[4,]

  return(list(vjind=vjind,v.alpha=v.alpha, v.beta=v.beta, a.mu=a.mu, a.s=a.s, b.mu=b.mu, b.s=b.s, k.mu=k.mu,k.s=k.s, v.eta=v.eta, v.etak=v.etak, v.z=v.z, v.muk=v.muk,jump.p=c(alpha.p, beta.p,eta.p,etak.p, z.p,muk.p)))
}


f.metro<-function(data, n.iter=3000, m.iter=3, n.thin=10,ls.dim=3,total.prop,z.pos.init,muk.fix){
  require(movMF)
  mc.data<-data
  data.n<-nrow(mc.data)
  print("data matrix dimensions are:")
  print(dim(mc.data))
  data.j<-ncol(mc.data)
  n.keep<-n.iter/n.thin
  sims<-array(NA, c(n.keep, m.iter, data.n+data.j+4+data.j+1))
  sims.latent<-array(NA,c(n.keep,m.iter,data.n+data.j,ls.dim))
  dimnames(sims)<-list(NULL, NULL, c(paste("alpha", 1:data.n, sep=''),
                                     paste("beta", 1:data.j, sep=''),
                                     'a.mu', 'a.s', 'b.mu', 'b.s',paste("etak", 1:data.j, sep=''),"eta"))
  dimnames(sims.latent)<-list(NULL, NULL, c(paste("respls", 1:data.n, sep=''),
                                            paste("subpopls", 1:data.j, sep='')),NULL)
  sims.p<-array(NA, c(n.keep, m.iter, data.n+data.j+data.n+data.j+data.j+1))
  dimnames(sims.p)<-list(NULL, NULL, c(paste("alpha", 1:data.n, sep=''),
                                       paste("beta", 1:data.j, sep=''),"eta",paste("etak", 1:data.j, sep=''),
                                       paste("resp", 1:data.n, sep=''),paste("subpop", 1:data.j, sep='')))
  for(m in 1:m.iter){
    ## initialization
    mc.v.alpha<-rnorm(data.n,log(50),2)
    mc.v.beta<-rnorm(data.j,log(.001),1)
    mc.v.beta[exp(mc.v.beta)>1]<-log(.03)
    mc.a.mu<-log(50)
    mc.a.s<-1
    mc.b.mu<-log(.001)
    mc.b.s<-1
    mc.v.z<-z.pos.init
    mc.v.muk<-matrix(rnorm(data.j*ls.dim),data.j,ls.dim)
    mc.v.muk<-mc.v.muk/sqrt(rowSums(mc.v.muk^2))
    mc.v.muk[muk.fix.ind[1],]<-muk.fix[1,]
    mc.v.muk[muk.fix.ind[2],]<-muk.fix[2,]
    mc.v.muk[muk.fix.ind[3],]<-muk.fix[3,]
    mc.v.muk[muk.fix.ind[4],]<-muk.fix[4,]
    
    mc.v.eta<-1
    mc.v.etak<-runif(data.j,0.001,50)
    mc.k.mu<-log(40)
    mc.k.s<-log(10)
    ad.jump<-rep(mc.a.s, data.n)
    bd.jump<-rep(mc.b.s, data.j)
    eta.jump=4
    etak.jump<-bd.jump
    muk.jump<-runif(length(bd.jump),150,300)
    z.jump<-runif(nrow(data),1,10)
    mc.j.ind<-rep(0,length(muk.jump))
    
    
    ### jump's are the variance of the jump distribution.
    last.50p<-array(NA, c((length(sims[1,1,])+length(sims.latent[1,1,,1])-4), 50))
    p.ct<-0
    for(t in 1:n.iter){
      temp<-f.update(mc.data, mc.v.alpha, mc.v.beta, mc.v.eta, mc.v.etak, mc.v.z, mc.v.muk, mc.a.mu, mc.a.s, mc.b.mu, mc.b.s, ad.jump, bd.jump,eta.jump,etak.jump,muk.jump,z.jump,mc.k.mu,mc.k.s,mc.j.ind,total.prop,muk.fix,ls.dim)
      mc.v.alpha<-temp$v.alpha
      mc.v.beta<-temp$v.beta
      mc.v.omega<-temp$v.omega
      mc.a.mu<-temp$a.mu
      mc.a.s<-temp$a.s
      mc.b.mu<-temp$b.mu
      mc.b.s<-temp$b.s
      mc.v.z<-temp$v.z
      mc.v.muk<-temp$v.muk
      mc.v.eta<-temp$v.eta
      mc.v.etak<-temp$v.etak
      mc.k.mu<-temp$k.mu
      mc.k.s<-temp$k.s
      mc.j.ind<-temp$vjind
      p.ct<-p.ct+1
      last.50p[,p.ct]<-temp$jump.p
      if(p.ct==50){
        print("summary of degree distribution")
        print(summary(exp(mc.v.alpha)))
        cat("a.s=",mc.a.s,"\n")
        jump.d<-c(ad.jump, bd.jump,eta.jump,etak.jump,z.jump,muk.jump)
        p.ct<-0
        p.mean<-rowMeans(last.50p, na.rm=T)
        print("jump prob respondent latent positions")
        print(summary(p.mean[(1+data.n+data.j+data.j+1):(1+data.n+data.j+data.j+data.n)]))
        print("etak jumping probability")
        print(round(p.mean[(1+data.n+data.j+1):(1+data.n+data.j+data.j)],3))
        print("current values of etak")
        print(mc.v.etak)
        jump.d.keep<-jump.d
        jump.d<-pmax(pmin(log(0.4)*jump.d/log(p.mean), 10*jump.d), 0.0001)

        ad.jump<-jump.d[1:data.n]
        bd.jump<-jump.d[(data.n+1):(data.n+data.j)]
        eta.jump<-jump.d[(1+data.n+data.j)]
        print("eta.jump prob")
        print(p.mean[(1+data.n+data.j)])
        print("current eta")
        print(mc.v.eta)

        etak.jump<-jump.d[(1+data.n+data.j+1):(1+data.n+data.j+data.j)]
        dtmp<-1/((log(.4)/log(p.mean))*(1/jump.d.keep))
        jump.d.keep<-pmax(pmin(1/((log(.4)/log(p.mean))*(1/jump.d.keep)), 750), 0.0001)
        muk.jump<-jump.d.keep[(1+data.n+data.j+data.j+data.n+1):(1+data.n+data.j+data.j+data.n+data.j)]
        z.jump<-jump.d.keep[(1+data.n+data.j+data.j+1):(1+data.n+data.j+data.j+data.n)]
      }
      if (t%%n.thin==0){
        print(c(t,m, date()))
        sims[t/n.thin,m,] <- c (mc.v.alpha, mc.v.beta, mc.a.mu, mc.a.s, mc.b.mu, mc.b.s,mc.v.etak,mc.v.eta)
        sims.latent[t/n.thin,m,,] <- rbind(mc.v.z,mc.v.muk)
        sims.p[t/n.thin,m,]<- temp$jump.p
      }
    }
  }
  return(list(sims=sims, sims.p=sims.p, sims.latent=sims.latent, data=mc.data, data.n=data.n, data.j=data.j, n.iter=n.iter, m.iter=m.iter, n.thin=n.thin))
}

