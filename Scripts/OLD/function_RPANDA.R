

full.phylo=tree
ana.events=ana_events_tables[[i]]
clado.events=clado_events_tables[[i]]
stratified=FALSE
map=group.map2
data=M
trim.class=subgroup
model="MC"
par=NULL
method="Nelder-Mead"
bounds=NULL


fit_t_comp_subgroup<-function(full.phylo,ana.events,clado.events,stratified=FALSE,map,data,trim.class,model=c("MC","DDexp","DDlin"),par=NULL,method="Nelder-Mead",bounds=NULL){
  
  if(is.null(names(data))){stop("data missing taxa names")}
  if(!is.null(dim(data))){stop("data needs to be a single trait")}
  is_tip <- full.phylo$edge[,2] <= length(full.phylo$tip.label)
  if(sum(diff(full.phylo$edge[is_tip, 2])<0)>0){ stop('fit_t_comp_subgroup cannot be used with ladderized full.phylogenies')}
  
  if(is.null(bounds[["lower"]]) & is.null(bounds[["upper"]])){
    bounds$lower = -Inf
    bounds$upper = Inf
  }
  
  GeoByClassObject<-CreateGeobyClassObject(full.phylo,map,trim.class,ana.events,clado.events,stratified=stratified)
  
  phylo<-GeoByClassObject$map
  
  root.trimmed.phylo<-max(nodeHeights(phylo))
  root.data<-max(nodeHeights(drop.tip.simmap(phylo,phylo$tip.label[which(!phylo$tip.label%in%names(data))])))
  
  if(round(root.trimmed.phylo,5)!=round(root.data,5)){stop("error where root of trimmed simmap and root of target clade don't match")}
  #if this above error is triggered, need to use scripts written by JPD using JC mvMORPH approach in summer 2018
  
  if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
  geo.object<-GeoByClassObject$geo.object
  #geo.sorted<-.resortGeoObject(phylo,geo.object) 
  
  if(length(geo.object$geography.object)<phylo$Nnode){stop("geography object cannot have more or fewer components than internode intervals in phylo")}
  if(length(data)>length(phylo$tip.label)){stop("error: some tips missing from pruned simmap")}
  if(!all(names(data)%in%phylo$tip.label)){stop("error: some tips missing from pruned simmap")}
  
  if(is.null(par)){par<-c(log(sqrt(var(data)/max(nodeHeights(extract.clade(phylo,getMRCA(phylo,names(data))))))),0)}
  
  if(model=="MC"){
    opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,data=data,model="MC",method=method, lower=bounds$lower, upper=bounds$upper)
    sig2 = exp(opt$par[1])^2
    S = -abs(opt$par[2])
    z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="MC",par=opt$par,return.z0=TRUE)
    results<-list(model = model, LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = sig2, S = S, z0 = as.numeric(z0), convergence = opt$convergence)
    return(results)
  }
  if(model=="DDexp"){
    opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,data=data,model="DDexp",method=method, lower=bounds$lower, upper=bounds$upper)
    sig2 = exp(opt$par[1])^2
    r = opt$par[2]
    z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="DDexp",par=opt$par,return.z0=TRUE)
    results<-list(model = model, LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = sig2, r = r, z0 = as.numeric(z0), convergence = opt$convergence)
    return(results)
  }
  if(model=="DDlin"){
    geography.matrix<-geo.object$geography.object
    maxN<-max(vapply(geography.matrix,function(x)max(rowSums(x)),1))
    opt<-optim(par,likelihood_subgroup_model,phylo=phylo,geography.object=geo.object,maxN=maxN,data=data,model="DDlin",method=method, lower=bounds$lower, upper=bounds$upper)
    sig2 = exp(opt$par[1])^2
    b = opt$par[2]
    z0=likelihood_subgroup_model(data=data,phylo=phylo,geography.object=geo.object,model="DDlin",par=opt$par,return.z0=TRUE,maxN=maxN)
    results<-list(model = model, LH = -opt$value, aic = (2*3 - 2*(-opt$value)), aicc = (2*3 - 2*(-opt$value))+((2*3*(3+1))/(length(phylo$tip.label)-3-1)), free.parameters = 3, sig2 = sig2, b = b, z0 = as.numeric(z0), convergence = opt$convergence)
    return(results)
  }
  
}
# 
# phylo=full.phylo
# simmap=map
# trim.class
# ana.events
# clado.events
# stratified=stratified


CreateGeobyClassObject<-function(phylo,simmap,trim.class,ana.events,clado.events,stratified=FALSE,rnd=5){
  
  trc=trim.class
  
  new.map<-.trimSimmap(simmap,trc)
  
  class.object<-CreateClassObject(new.map)
  geo.object<-CreateGeoObject_BioGeoBEARS(full.phylo=phylo,trimmed.phylo=new.map,ana.events,clado.events,stratified=stratified)
  
  geot<-round(geo.object$times,rnd)
  clat<-round(class.object$times,rnd)
  nodeDist<-sort(unique(c(geot,clat)))
  nodeDiff<-diff(nodeDist)
  if(any(nodeDiff<= (2*(10^-rnd)))){stop("potential rounding error, two time bins very similar, try changing rnd digits")}
  
  nat<-list()
  
  u<-0
  y<-0
  
  for(i in 1:length(nodeDiff)){
    
    if((nodeDist[i]%in%geot) && (nodeDist[i]%in%clat)){ #if timing is the same for both
      u = u+1
      y = y+1
      tr.vec<-rep(0,nrow(class.object$class.object[[u]]))
      tr.vec[which(class.object$class.object[[u]][,2]%in%trc)]<-1
      hold.mat<-geo.object$geography.object[[y]]*tr.vec%*%t(tr.vec)
      diag(hold.mat)<-1
      nat[[i]]<-hold.mat
    }
    if((nodeDist[i]%in%geot) && (!nodeDist[i]%in%clat)){ #this means that geo.object changes but class object doesn't
      y = y+1
      tr.vec<-rep(0,nrow(class.object$class.object[[u]]))
      tr.vec[which(class.object$class.object[[u]][,2]%in%trc)]<-1
      hold.mat<-geo.object$geography.object[[y]]*tr.vec%*%t(tr.vec)
      diag(hold.mat)<-1
      nat[[i]]<-hold.mat		
    }
    if((!nodeDist[i]%in%geot) && (nodeDist[i]%in%clat)){ #this means that geo.object changes but class object doesn't
      u = u+1
      tr.vec<-rep(0,nrow(class.object$class.object[[u]]))
      tr.vec[which(class.object$class.object[[u]][,2]%in%trc)]<-1
      hold.mat<-geo.object$geography.object[[y]]*tr.vec%*%t(tr.vec)
      diag(hold.mat)<-1
      nat[[i]]<-hold.mat		
    }
  }
  
  coll.vector<-vector()
  count=1
  geography.matrix<-list()
  for(i in 1:length(nat)){
    int.mat<-nat[[i]] 
    if(count>1 && !is.character(all.equal(nat[[i]],geography.matrix[[count-1]]))){coll.vector<-c(coll.vector,i)} else{	
      geography.matrix[[count]]<-int.mat
      count=count+1
    }
  }
  nodeDist<-nodeDist[which(!nodeDist%in%nodeDist[coll.vector])]
  nodeDist<-nodeDist-min(nodeDist)
  nodeDiff<-diff(nodeDist)	
  
  return(list(map=new.map,geo.object=list(geography.object=geography.matrix,times=nodeDist,spans=nodeDiff)))#new phylo object, #new times, #new spans, #new geo object
}


.trimSimmap<-function(map,trim.class){
  trc=trim.class
  smap<-map
  
  repeat{
    
    #print("start")
    
    tot.len<-length(smap$tip.label)
    node = tot.len + 1
    tbdropped<-c()
    
    while(node <= (tot.len+smap$Nnode)){
      descL<-smap$edge[which(smap$edge[,1]==node),2][1]
      Lclade<-unique(c(which(smap$edge[,1]==node)[1],which(smap$edge[,2]%in%getDescendants(smap,descL))))
      
      if(!any(sapply(smap$maps[Lclade],function(x)any(names(x)%in%trc)))){ #left descendants DO NOT contain target trait
        
        if(length(Lclade)==1){ 	#if only descendant is tip, drop entire branch
          tbdropped<-c(tbdropped,smap$tip.label[descL])					
        }
        
      }
      
      node = node + 1
      
    }	
    
    smap<-drop.tip.simmap(smap,tbdropped)
    smapL<-smap
    
    #print(node)
    
    tot.len<-length(smap$tip.label)
    node = tot.len + 1
    tbdropped<-c()
    
    while(node <= (tot.len+smap$Nnode)){
      descR<-smap$edge[which(smap$edge[,1]==node),2][2]
      Rclade<-unique(c(which(smap$edge[,1]==node)[2],which(smap$edge[,2]%in%getDescendants(smap,descR))))
      
      if(!any(sapply(smap$maps[Rclade],function(x)any(names(x)%in%trc)))){ #right descendants DO NOT contain target trait
        
        if(length(Rclade)==1){ 	#if only descendant is tip, drop entire branch
          tbdropped<-c(tbdropped,smap$tip.label[descR])
          
        } 
        
      }
      
      node = node +1
      
    }
    
    smap<-drop.tip.simmap(smap,tbdropped)
    
    #print("error2")
    
    if(length(smapL$tip.label)==length(smap$tip.label)){
      break
    }
  }
  
  return(smap)
}	

CreateClassObject<-function(simmap,rnd=5){
  if(any(grepl("___",simmap$tip.label))){stop("script will not work with '___' in tip labels; remove extra underscores")}
  paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
  paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
  nodeDist<-vector(mode = "numeric", length = simmap$Nnode)
  totlen<-length(simmap$tip.label)
  root <-totlen  + 1
  heights<-nodeHeights(simmap)
  for (i in 1:dim(simmap$edge)[1]){
    nodeDist[[simmap$edge[i, 1] - totlen]] <- heights[i]
  }
  nodeDist<-c(nodeDist,max(heights))
  nodeDiff<-diff(nodeDist)
  flag=0
  if(sum(nodeDiff<0)>0){  ##this loop renumbers the nodes if trees nodes are not placed in sequential order
    node.order<-match(rank(heights[,1],ties.method="min"),seq(1, by = 2, len = simmap$Nnode))
    node.order<-node.order+totlen
    old.edge<-simmap$edge
    old.map<-simmap
    simmap$edge[,1]<-node.order
    for(j in 1:length(simmap$edge[,2])){
      if(simmap$edge[j,2]>totlen){
        #match number order in old edge
        #lookup value in new edge
        #replace with value
        simmap$edge[j,2]<-simmap$edge[,1][match(simmap$edge[j,2],old.edge[,1])]
      }
    }
    nodeDist<-vector()
    for (i in 1:dim(simmap$edge)[1]){
      nodeDist[[simmap$edge[i, 1] - totlen]] <- heights[i]
    }
    nodeDist<-c(nodeDist,max(heights))
    nodeDiff<-diff(nodeDist)
    flag=1
  }
  mat<-matrix(nrow=0, ncol=3)
  counter_three_letters <- 0
  for(i in 1:simmap$Nnode){
    other<-simmap$edge[simmap$edge[,1]==i+totlen, 2]
    for(b in other){
      int<-matrix(ncol=3)
      int[1]<-i+totlen
      if(b>totlen){
        counter_three_letters <- counter_three_letters + 1
        int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
        int[3]<-b
      } else {
        int[2]<-simmap$tip.label[[b]]
        int[3]<-0 
      }
      mat<-rbind(mat,int)
    }
  }
  if(any(class(simmap)=="phylo")){
    for(i in 1:length(mat[,1])){
      if(mat[i,3]==0){
        mat[i,3]<-as.character(match(mat[i,2],simmap$tip.label))
      }}	
    maps.object<-simmap$maps
    brchange<-vapply(maps.object,function(x) any(length(names(x))>1),1)
    if(sum(brchange)!=0){
      newtimes<-vector()
      for(i in 1:length(brchange)){
        if(brchange[i]!=0){
          #look up branch in stochastic map
          intlen<-length(maps.object[[i]])
          inttime<-vector()
          for(j in 1:(intlen-1)){
            inttime<-c(inttime,maps.object[[i]][j])
            nt<-as.numeric(heights[i,1]+sum(inttime))
            newtimes<-c(newtimes,nt)
          }
        }
      }
      old.Dist<-nodeDist
      old.Diff<-nodeDiff
      nodeDist<-sort(c(nodeDist,newtimes))
      nodeDiff<-diff(nodeDist)
    }
    if(any(nodeDiff<= (2*(10^-rnd)))){stop("CreateClassObject.R: potential rounding error, two time bins very similar")}
    nat<-list()
    nodecount=1
    for(i in 1:length(nodeDiff)){
      if(i==1){
        hold.m<-mat[as.numeric(mat[,1])<=(totlen+nodecount),c(2,3)]
        int<-dim(hold.m)[1]
        for(m in 1:int){
          hold.m[m,2]<-names(maps.object[[which(simmap$edge[,2]==as.numeric(hold.m[m,2]))]])[1]
        }
        nat[[i]]<-hold.m
      } else {
        if(nodeDist[i]%in%old.Dist){ #a new species appears, this is a node
          P<-mat[as.numeric(mat[,1])<=(totlen+nodecount),c(2,3)]
          hold.m<-rbind(P[as.numeric(P[,2])<=totlen,],P[as.numeric(P[,2])>(totlen+nodecount),])
          int<-dim(hold.m)[1]
          for(m in 1:int){
            iden<-which(simmap$edge[,2]==as.numeric(hold.m[m,2]))
            if(brchange[iden]==0 || !(hold.m[m,1]%in%nat[[i-1]][,1])){
              hold.m[m,2]<-names(maps.object[[iden]])[1]
            } else{
              num=1
              while(round(nodeDist[i+1]-old.Dist[simmap$edge[iden,1]-totlen],6)>round(sum(maps.object[[iden]][1:num]),6)){
                num=num+1}
              hold.m[m,2]<-names(maps.object[[iden]])[num]
            }
          }
          nat[[i]]<-hold.m
        } else{ # this is a branch change
          P<-mat[as.numeric(mat[,1])<=(totlen+nodecount),c(2,3)]
          hold.m<-rbind(P[as.numeric(P[,2])<=totlen,],P[as.numeric(P[,2])>(totlen+nodecount),])
          int<-dim(hold.m)[1]
          for(m in 1:int){
            iden<-which(simmap$edge[,2]==as.numeric(hold.m[m,2]))
            if(brchange[iden]==0 | !(hold.m[m,1]%in%nat[[i-1]][,1])){
              hold.m[m,2]<-names(maps.object[[iden]])[1]
            } else{					
              num=1
              while((round(nodeDist[i+1]-old.Dist[simmap$edge[iden,1]-totlen],rnd)-round(sum(simmaps.object[[iden]][1:num]),rnd))>(2*(10^-rnd))){
                num=num+1}
              hold.m[m,2]<-names(simmaps.object[[iden]])[num]
            }
          }
          nat[[i]]<-hold.m
        }
      }
      if(nodeDist[i+1]%in%old.Dist){nodecount=nodecount+1}
    }	
  }
  return(list(class.object=nat,times=round(nodeDist,8),spans=nodeDiff))
  
}

CreateGeoObject_BioGeoBEARS<-function(full.phylo,trimmed.phylo=NULL,ana.events,clado.events,stratified=FALSE){
  
  if(stratified){
    if(is.null(trimmed.phylo)){
      clado_events_tables<-list()
      clado_events_tables[[1]]<-clado.events
      smap<-.stratified_BGB_to_tables(full.phylo,clado_events_tables,1)
      x<-.CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=full.phylo,ana.events=smap$ana.int,clado.events=smap$clado.int,nat.only=FALSE)
      return(x)
    }else{
      clado_events_tables<-list()
      clado_events_tables[[1]]<-clado.events
      smap<-.stratified_BGB_to_tables(full.phylo,clado_events_tables,1)
      x<-.CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=trimmed.phylo,ana.events=smap$ana.int,clado.events=smap$clado.int,nat.only=FALSE)
      return(x)		
    }		
  } else{
    if(is.null(trimmed.phylo)){
      x<-.CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=full.phylo,ana.events=ana.events,clado.events=clado.events,nat.only=FALSE)
      return(x)
    } else{
      x<-.CreateBioGeoB_Object_subclade(anc.phylo=full.phylo,subclade.phylo=trimmed.phylo,ana.events=ana.events,clado.events=clado.events,nat.only=FALSE)
      return(x)
    }
  }
}

.CreateBioGeoB_Object_subclade<-function(anc.phylo,subclade.phylo,ana.events,clado.events,nat.only=FALSE){
  phylo<-subclade.phylo
  subclade.tips<-phylo$tip.label
  paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
  paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS
  totlen<-length(phylo$tip.label)
  root <-totlen  + 1
  heights<-nodeHeights(phylo)
  totheight=max(nodeHeights(phylo))
  totheight.anc=max(nodeHeights(anc.phylo))
  subclade.root=getMRCA(anc.phylo,subclade.tips)
  tips = clado.events$node[which(clado.events$label%in%subclade.tips)]
  nn<-unique(apply(combn(tips,2),2,function(x)getMRCA(anc.phylo,x)))
  M<-mrca(anc.phylo,full=TRUE)
  desc=unique(as.vector(M[which(rownames(M)%in%c(tips,nn)),]))
  desc=desc[c(which(desc>=subclade.root),which(desc<=length(anc.phylo$tip.label)))]
  subnodes=c(desc[which(desc>length(anc.phylo$tip.label))])
  ghost.nodes=desc[which(!desc%in%c(nn,tips))]
  ana.events<-ana.events[which(ana.events$node%in%desc[which(desc!=subclade.root)]),]
  clado.events<-clado.events[which(clado.events$node%in%subnodes),]
  nodes=unique(round(clado.events[,9],8))
  events=unique(round(totheight.anc-ana.events$abs_event_time,8))
  nodeDist=sort(c(nodes,events,totheight.anc))
  nodeDiff=diff(nodeDist)
  old.phylo<-phylo
  old.edge<-anc.phylo$edge[which(anc.phylo$edge[,1]%in%nn),]
  other.edge<-phylo$edge
  if(!is.null(phylo$node.label)){phylo$node.label<-NULL}
  old.labels<-as.numeric(names(sort(branching.times(phylo),decreasing=TRUE)))
  if(any(diff(old.labels)!=1)){ #if nodes are not in sequential order, this renames them so that they are
    checkmat<-cbind(old.labels,seq(root,length(phylo$tip.label)+phylo$Nnode))
    for(j in 1:phylo$Nnode){phylo$edge[which(other.edge==checkmat[j,1])]<-checkmat[j,2]}
  }	
  mat<-matrix(nrow=0, ncol=3)
  counter_three_letters <- 0
  for(i in 1:phylo$Nnode){
    other<-phylo$edge[phylo$edge[,1]==i+totlen, 2]
    for(b in other){
      int<-matrix(ncol=3)
      int[1]<-i+totlen
      if(b>totlen){
        counter_three_letters <- counter_three_letters + 1
        int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
        int[3]<-b
      } else {
        int[2]<-phylo$tip.label[[b]]
        int[3]<-0 
      }
      mat<-rbind(mat,int)
    }
  }
  desc.mat<-matrix(nrow=length(desc), ncol=2)
  desc.mat[,1]<-desc
  for(i in 1:length(desc)){
    if(desc.mat[i,1]%in%ghost.nodes){
      flag=0
      j=1
      while(flag==0){
        term=anc.phylo$edge[which(anc.phylo$edge[,1]==desc.mat[i,1]),2][j]
        if(term%in%desc){
          flag=1
        }else{
          j=2
        }}
      if(term%in%tips){
        desc.mat[i,2]<-which(anc.phylo$tip.label[term]==phylo$tip.label)
      } 
      if(term%in%nn){
        desc.mat[i,2]<-phylo$edge[which(old.edge[,1]==term),1][1]
      }
      term2=term 
      while(term2%in%ghost.nodes){
        flag=0
        j=1
        while(flag==0){
          term=anc.phylo$edge[which(anc.phylo$edge[,1]==term2),2][j]
          if(term%in%desc){
            flag=1
          }else{
            j=2
          }}
        if(term%in%tips){
          desc.mat[i,2]<-which(anc.phylo$tip.label[term]==phylo$tip.label)
        } 
        if(term%in%nn){
          desc.mat[i,2]<-phylo$edge[which(old.edge[,1]==term),1][1]
        } 
        term2=term
      }
    } 
    if(desc.mat[i,1]%in%tips){
      desc.mat[i,2]<-which(anc.phylo$tip.label[desc.mat[i,1]]==phylo$tip.label)
    }
    if(desc.mat[i,1]%in%nn){
      desc.mat[i,2]<-phylo$edge[which(old.edge[,1]==desc.mat[i,1]),1][1]
    }
  }
  
  
  #while(any(desc.mat[,2]%in%ghost.nodes)){
  #	nxt<-desc.mat[which(desc.mat[,2]%in%ghost.nodes),2]
  #	desc.mat[which(desc.mat[,2]%in%ghost.nodes),2]<-desc.mat[match(nxt ,desc.mat[,1]),2]
  #}
  
  nat<-list()
  nodecount=1
  for(i in 1:length(nodeDiff)){
    if(nodeDist[i]%in%nodes){ #a new species appears, this is a node
      if(clado.events[which(round(clado.events[,9],8)==round(nodeDist[i],8)),1]%in%nn){## node has two descendants because of subtree trimming
        IN<-vector()
        node<-totlen+nodecount
        P<-mat[as.numeric(mat[,1])<=(node),c(2,3)]
        IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(node),1])
        #need to find a way to look up old node
        #oldnode=old.edge[which(phylo$edge[,1]==node)][1]
        delbr<-clado.events[which(round(clado.events[,9],8)==round(nodeDist[i],8)),1]
        left=clado.events[which(clado.events[,1]==delbr),15]
        left=desc.mat[which(desc.mat[,1]==left),2]
        right=clado.events[which(clado.events[,1]==delbr),16]
        right=desc.mat[which(desc.mat[,1]==right),2]
        m<-regexpr("->",clado.events[which(clado.events[,1]==delbr),20])
        regs<-regmatches(clado.events[which(clado.events[,1]==delbr),20],m,invert=TRUE)[[1]][2]
        p<-regexpr(",",regs)	
        sides<-regmatches(regs,p,invert=TRUE)
        lside<-sides[[1]][1]
        rside<-sides[[1]][2]
        if(i==1){
          geo.vector<-vector(length=length(IN))
          if(left<=totlen){
            geo.vector[which(IN==phylo$tip.label[left])]<-lside
          }else{
            geo.vector[which(IN==mat[,2][mat[,3]==left])]<-lside
          }
          if(right<=totlen){
            geo.vector[which(IN==phylo$tip.label[right])]<-rside
          }else{
            geo.vector[which(IN==mat[,2][mat[,3]==right])]<-rside
          }
          nat[[i]]<-cbind(IN,geo.vector)
        } else {
          ogv<-geo.vector			##need to update old geo.vector
          geo.vector<-vector(length=length(IN))
          terms<-which(nat[[i-1]][,1]%in%IN) #this should give the item numbers of old things
          geo.vector[match(nat[[i-1]][,1][terms],IN)]<-ogv[terms]
          if(left<=totlen){
            geo.vector[which(IN==phylo$tip.label[left])]<-lside
          }else{
            geo.vector[which(IN==mat[,2][mat[,3]==left])]<-lside
          }
          if(right<=totlen){
            geo.vector[which(IN==phylo$tip.label[right])]<-rside
          }else{
            geo.vector[which(IN==mat[,2][mat[,3]==right])]<-rside
          }
          nat[[i]]<-cbind(IN,geo.vector)			
        }} else{
          IN<-nat[[i-1]][,1]
          delbr<-clado.events[which(round(clado.events[,9],8)==round(nodeDist[i],8)),1]
          left=clado.events[which(clado.events[,1]==delbr),15]
          right=clado.events[which(clado.events[,1]==delbr),16]
          m<-regexpr("->",clado.events[which(clado.events[,1]==delbr),20])
          regs<-regmatches(clado.events[which(clado.events[,1]==delbr),20],m,invert=TRUE)[[1]][2]
          p<-regexpr(",",regs)	
          sides<-regmatches(regs,p,invert=TRUE)
          lside<-sides[[1]][1]
          rside<-sides[[1]][2]
          if(left%in%desc){ #left branch should be recorded
            hold<-desc.mat[which(desc.mat[,1]==left),2]
            if(hold<=totlen){
              elno<-which(IN==phylo$tip.label[hold])
            }else{
              elno<-which(IN==mat[,2][mat[,3]==hold])
            }
            geo.vector[elno]<-lside
          } else{#right branch should be recorded
            hold<-desc.mat[which(desc.mat[,1]==right),2]
            if(hold<=totlen){
              elno<-which(IN==phylo$tip.label[hold])
            }else{
              elno<-which(IN==mat[,2][mat[,3]==hold])
            }
            geo.vector[elno]<-rside
          }
          nat[[i]]<-cbind(IN,geo.vector)			
        }} else{ # this is a branch change
          IN<-nat[[i-1]][,1]
          #identify which branch changes
          ##note should eventually move this up so subtraction doesn't happen every time(for speed up)
          delbr<-which(nodeDist[i]==round(totheight.anc-ana.events$abs_event_time,8))
          if(length(delbr)>1){stop("more than one events at same time")}
          hold<-desc.mat[which(desc.mat[,1]==ana.events[delbr,]$nodenum_at_top_of_branch),2]
          if(hold<=totlen){
            elno<-which(IN==phylo$tip.label[hold])
          }else{
            elno<-which(IN==mat[,2][mat[,3]==hold])
          }
          #identify new state
          geo.vector[elno]<-ana.events[delbr,]$new_rangetxt
          nat[[i]]<-cbind(IN,geo.vector)			
        }
    if(nodeDist[i+1]%in%nodes && clado.events[which(round(clado.events[,9],8)==round(nodeDist[i+1],8)),1]%in%nn){nodecount=nodecount+1}
    if(any(nat[[i]][,2]==FALSE)){stop("ERROR:'FALSE' recorded as range")}
  }			
  
  if(nat.only==TRUE){
    return(nat[[length(nat)]])
  } else{
    
    coll.vector<-vector()
    count=1
    geography.matrix<-list()
    for(i in 1:length(nodeDiff)){ 
      timei<-unlist(nat[[i]])
      var.list<-timei[,1]
      geo.list<-strsplit(timei[,2],"")
      len=length(var.list)
      int.mat<-matrix(nrow=len,ncol=len)
      rownames(int.mat)<-var.list
      colnames(int.mat)<-var.list
      diag(int.mat)<-1
      for(j in 1:length(var.list)){
        for(k in 1:length(var.list)){
          if(j>k){
            int.mat[j,k]<-ifelse(any(geo.list[[j]]%in%geo.list[[k]]),1,0)
          }
        }
      }
      int.mat[upper.tri(int.mat)==TRUE]<-t(int.mat)[upper.tri(t(int.mat))==TRUE]
      if(count>1 && !is.character(all.equal(int.mat,geography.matrix[[count-1]]))){coll.vector<-c(coll.vector,i)} else{	
        geography.matrix[[count]]<-int.mat
        count=count+1
      }
    }
    nodeDist<-nodeDist[which(!nodeDist%in%nodeDist[coll.vector])]
    nodeDist<-nodeDist-min(nodeDist)
    nodeDiff<-diff(nodeDist)	
    return(list(geography.object=geography.matrix,times=nodeDist,spans=nodeDiff))
  }
}
