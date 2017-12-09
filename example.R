library(igraph)
load('/Users/Mengjie/Desktop/social_network/data/zdata.RData')
load('/Users/Mengjie/Desktop/social_network/data/znets.RData')
load('/Users/Mengjie/Desktop/social_network/data/hh.features.new.RData')
load("/Users/Mengjie/Desktop/social_network/code_share/distance.all_2.RData")

network.list=list()
degrees.list=list()
centrality.list=list()
closeness.list=list()
betweenness.list=list()
count.triangle.list=list()
support.list=list()
max.eigenvalue.list=NULL
average.path.length.list=NULL
diameter.list=NULL
fraction.giant.component=NULL
num.components.list=NULL
for (i in c(1:12,14:21,23:77)){
  reluse=1
  villagei.adj=as.matrix(znets[[1]][[i]][[1]][[reluse]][[1]])  
  for (reluse in c(2:11,13)){
    villagei.adj=pmax(villagei.adj,as.matrix(znets[[1]][[i]][[1]][[reluse]][[1]]),na.rm=T)
  }
  n=dim(villagei.adj)[1]
  network.list=c(network.list,list(villagei.adj))
  graphi=graph.adjacency(villagei.adj,mode='undirected')
  degrees.list=c(degrees.list,list(rowSums( villagei.adj)))
  centrality=evcent(graphi,scale=F)$vector
  centrality.list=c(centrality.list,list(centrality))
  closeness.list=c(closeness.list,list(closeness(graphi)))
  betweenness.list=c(betweenness.list,list(betweenness(graphi)))
  count.triangle.list=c(count.triangle.list,list(count_triangles(graphi)))
  max.eigenvalue.list=c(max.eigenvalue.list,evcent(graphi,scale=F)$value)
  average.path.length.list=c(average.path.length.list,mean_distance(graphi))
  diameter.list=c(diameter.list,diameter(graphi))
  fraction.giant.component=c(fraction.giant.component,max(components(graphi)$csize)/n)
  num.components.list=c(num.components.list,count_components(graphi))
  support.t=NULL
  for (ind in 1:n){
    nb=which(graphi[ind,]==1)
    if(length(nb)>=2){support.t=c(support.t,length(which(rowSums(graphi[nb,nb])>0)))}
    else{support.t=c(support.t,0)}
  }
  support.list=c(support.list,list(support.t))
}

zdata.one=zdata[!duplicated(zdata$calc_hh_id),]
zdata.one$village.id=as.numeric(c(substr(zdata.one$calc_hh_id[1:529],1,1),substr(zdata.one$calc_hh_id[530:dim(zdata.one)[1]],1,2)))
for(i in 1:length(zdata.one$village.id)){
  if(zdata.one$village.id[i]>22){zdata.one$village.id[i]=zdata.one$village.id[i]-2}
  else if(zdata.one$village.id[i]>13){zdata.one$village.id[i]=zdata.one$village.id[i]-1}
}

zdata.hh.id=NULL
for (vlg in 1:75){
  village=zdata.one[which(zdata.one$village.id==vlg),]
  village$hh.id=rep(NA,dim(village)[1])
  for (i in 1:dim(village)[1]){
    village$hh.id[i]=which(hh.features.list[[vlg]]$newhhid==village$calc_hh_id[i])
  }
  zdata.hh.id=c(zdata.hh.id,village$hh.id)
}
zdata.one$hh.ind=zdata.hh.id

data.use=subset(zdata.one,select=c(village.id,hh.ind,calc_hh_id,calc_gender,s2_2_num_twin_birth,s2_4_num_kas,s2_6_num_child_blw18,s2_8_num_10th,s2_10_num_ug_vlg,s2_12_num_4whelr_vlg,s2_14_vlg_tractor,s2_16_vlg_vhcle,s2_18_vlg_disease,s2_20_greatr_1wife_vlg,s2_22vlg_widow,s2_24vlg_leader))
data.use=data.use[-which(is.na(data.use),arr.ind=T)[,1],]
data.char=subset(zdata.one,select=c(village.id,hh.ind,calc_hh_id,calc_gender,s2_1_twin_birth,s2_3_kas,s2_5_child_blw18,s2_7_10th,s2_9_hh_ug,s2_11_hh_4whelr,s2_13_hh_tractor,s2_15_hh_vhcle,s2_17_hh_disease,s2_19_greatr_1wife_hh,s2_21hh_widow,s2_23hh_leader))
data.char=data.char[-which(is.na(data.char),arr.ind=T)[,1],]

network.list.subset=list()
degrees.list.subset=list()
centrality.list.subset=list()
closeness.list.subset=list()
betweenness.list.subset=list()
max.eigenvalue.list.subset=NULL
avg.dist.list.subset=list()
for (i in 1:75){
  network.list.subset=c(network.list.subset,list(network.list[[i]][data.use[which(data.use$village.id==i),]$hh.ind,data.use[which(data.use$village.id==i),]$hh.ind]))
  graphi=graph.adjacency(network.list.subset[[i]],mode='undirected')
  degrees.list.subset=c(degrees.list.subset,list(rowSums( network.list.subset[[i]])))
  centrality=evcent(graphi,scale=F)$vector
  centrality.list.subset=c(centrality.list.subset,list(centrality))
  max.eigenvalue.list.subset=c(max.eigenvalue.list.subset,evcent(graphi,scale=F)$value)
  closeness.list.subset=c(closeness.list.subset,list(closeness(graphi)))
  betweenness.list.subset=c(betweenness.list.subset,list(betweenness(graphi)))
  dist.table=as.matrix(distances(graphi))
  avg.dist=mean(1/dist.table[upper.tri(dist.table, diag = FALSE)])
  avg.dist.list.subset=c(avg.dist.list.subset,avg.dist)
}

data.use.constructARD=data.use
for(i in 1:75){
  data.char.vlg=data.char[which(data.char$village.id==i),]
  for(k in 1:12){
    if(length(which(data.char.vlg[,(k+4)]==1))>1){  
      data.use.constructARD[which(data.use.constructARD$village.id==i),k+4]=rowSums(network.list.subset[[i]][,which(data.char.vlg[,(k+4)]==1)])}
    else if(length(which(data.char.vlg[,(k+4)]==1))==1){
      data.use.constructARD[which(data.use.constructARD$village.id==i),k+4]=(network.list.subset[[i]][,which(data.char.vlg[,(k+4)]==1)])}
    else{
      data.use.constructARD[which(data.use.constructARD$village.id==i),k+4]=numeric(length(which(data.use.constructARD$village.id==i)))  
    }
  }
}

total.prop=NULL
x.axis=NULL
for (vlg in 1:75){
  villagei=data.char[which(data.char$village.id==vlg),]
  villagei[which(villagei<0,arr.ind=T)]=NA
  n=dim(villagei)[1]
  temp=sum(x.axis)
  for (k in c(4:10,12)){
    x.axis=c(x.axis,sum(as.numeric(villagei[,k+4]==1),na.rm = T)/length(!is.na(villagei[,k+4])))
  }
  total.prop=c(total.prop,sum(x.axis)-temp)
}

source('main.R')
g.sims=list()
for (vlg in 1:75){
  y=data.use.constructARD[data.use.constructARD$village.id==vlg,c(4:10,12)+4]
  y[which(y<0,arr.ind=T)]=NA
  y=as.matrix(y)
  muk.fix.ind=sample(1:8,size=4,replace=F)
  muk.fix=matrix(rnorm(12),nrow=4,ncol=3)
  muk.fix=sweep(muk.fix,MARGIN=1,1/sqrt(rowSums(muk.fix^2)),`*`)
  result=main(y=y,total.prop=total.prop[vlg],muk.fix=muk.fix,n.iter=3000, m.iter=3, n.thin=10,
              is.sample=TRUE,distance.matrix=distance.all[[vlg]],Knn.K=5)
  g.sims=c(g.sims,list(result))
  save(g.sims,file="g.sims_codeshare.RData")
}

#check random initial points

z=generateRandomInitial(500,3)

theta.fcn<-function(m){
  x<-m[1]
  y<-m[2]
  z<-m[3]
  theta<-asin(z/sqrt(x^2+y^2+z^2))
  phi<-atan2(x,y)
  return(list(theta=theta,phi=phi))
}

plot(1,1,xlim=c(-pi,pi),ylim=c(-pi/2,pi/2),type='n',axes=F,xlab='',ylab='',main='')
axis(side=1,at=c(-(.25+pi),-(.1+(pi/2)),0,(pi/2),pi),labels=c(0,expression(pi/2),expression(pi),expression(3*pi/2),expression(2*pi)),cex.axis=1.5,las=1)
axis(side=2,at=c(-(.12+pi/2),(-.05+(-pi/4)),0,pi/4,pi/2),labels=c(0,expression(pi/4),expression(pi/2),expression(3*pi/4),expression(pi)),cex.axis=1.5,las=1)
for(s in c(1:500)){
  zk.tmp<-theta.fcn(as.numeric(z[s,]))
  points(y=as.numeric(zk.tmp[1]),x=zk.tmp[2],col='red',pch=16)}


